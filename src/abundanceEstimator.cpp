/*
 * Copyright (C) 2018 Annika Seidel <annika.seidel@mpibpc.mpg.de>
 *
 */

#include <cstdio>
#include <fcntl.h>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>
#include <gatb/gatb_core.hpp>

#include "ToolInfo.h"
#include "options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "Lookuptable.h"
#include "CountProfile.h"
#include "preprocessing.h"
#include "processSequences.h"
#include "filehandling.h"

std::mutex sampleList_mutex;
std::chrono::milliseconds interval(100);



int writePopulationCoverage(CountProfile &countprofile, FILE *resultFile)
{
  /* estimate population coverage */
  float popCoverage = countprofile.calcPopulationCoverage() * countprofile.getCorrFactor();
  fprintf(resultFile, "%s\t%f\n", countprofile.getReadName(), popCoverage );
}

void thread_runner(int id, vector<string> *sampleList, std::string outdir,
                   Lookuptable *lookuptable, KmerTranslator *translator)
{
  std::string sampleFileName;
  string tempResultFileName;
  while(!sampleList->empty())
  {
    {
      while(!sampleList_mutex.try_lock() && !sampleList->empty())
      {
        std::this_thread::sleep_for(interval);
      }
      if(sampleList->empty())
      {
        fprintf(stderr, "Thread %d exiting, nothing more to do\n", id);
        return;
      }
      sampleFileName = sampleList->back();
      sampleList->pop_back();
      sampleList_mutex.unlock();

      tempResultFileName = outdir+basename(sampleFileName.c_str());

      std::cerr << "thread " << id << " working on "
                << sampleFileName << " outfilename: " << tempResultFileName
                << std::endl << std::flush;
    }
    processSeqFile(sampleFileName, tempResultFileName,
                    *lookuptable, *translator,writePopulationCoverage);

  }

}

void concatenate_read_files(vector<string> *sampleList,
                            std::string resultFileName,
                            std::string tmpOutDir)
{
  //Open Result File
  fprintf(stderr, "start concatenate read files to final results\n");
  char buf[BUFSIZ];
  size_t size;
  int resultFile = open(resultFileName.c_str(),
                  O_WRONLY | O_CREAT | O_TRUNC, 0660);
  if (resultFile == -1)
  {
    fprintf(stderr, "Error when trying to open file %s:\n%s\n",
            resultFileName.c_str(), strerror(errno));
    exit(EXIT_FAILURE);
  }
  /*
  // Write header
  if (write(resultFile, "#read_index\tpop_coverage\n", 25) == -1)
  {
    fprintf(stderr, "Error when writing to file %s:\n%s\n",
            resultFileName.c_str(), strerror(errno));
    close(resultFile);
    exit(EXIT_FAILURE);
  }*/

  //Iterate tmp resultFiles and copy to resultFile
  for(string sampleFileName: *sampleList)
  {
    string s_filename = tmpOutDir + basename(sampleFileName.c_str());

    int source = open(s_filename.c_str(), O_RDONLY, 0);
    if (source == -1)
    {
      fprintf(stderr, "ERROR: Opening result file %s failed:\n%s\n",
              s_filename.c_str(), strerror(errno));
      exit(EXIT_FAILURE);
    }
    while ((size = read(source, buf, BUFSIZ)) > 0)
    {
      if(write(resultFile, buf, size) == -1)
      {
        fprintf(stderr, "Error when writing to file %s: \n%s\n",
                resultFileName.c_str(), strerror(errno));
        close(source);
        close(resultFile);
        exit(EXIT_FAILURE);
      }
    }
    close(source);
    remove(s_filename.c_str());
  }
  close(resultFile);

}

void process_sampleList_threads(vector<string> *sampleList,
                                std::string resultFileName,
                                Lookuptable *lookuptable,
                                KmerTranslator *translator,
                                int num_threads)
{
  vector<string> *cp_sampleList = new vector<string>(*sampleList);
  //std::cout << "processing with " << num_threads << " threads" << std::endl
  //          << std::flush;
 
  //TODO: next lines as function! 
  string filename = resultFileName.c_str();
  size_t lastdot = filename.find_last_of(".");
    if (lastdot != std::string::npos)
      filename=filename.substr(0, lastdot);

  string outdir = string("tmp/tmp_") + basename(filename.c_str()) + "/";
  _mkdir(outdir);
  std::thread **threads = new std::thread*[num_threads];
  for(int i = 0; i < num_threads; i++)
  {
    threads[i] = new std::thread(thread_runner, i, cp_sampleList, outdir,
                                 lookuptable, translator);
  }
  for(int i = 0; i < num_threads; i++)
  {
    threads[i]->join();
    delete threads[i];
  }
  delete[] threads;
  concatenate_read_files(sampleList,resultFileName, outdir);
}



int abundanceEstimator(int argc, const char **argv, const ToolInfo* tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);
  printf("threads %u\n", opt.threads);
  printf("seqFile: %s\n", opt.seqFile.c_str());
  printf("kcFile: %s\n", opt.kcFile.c_str());
  printf("kmerWeight: %u\n", opt.kmerWeight);

  // TODO:check parameter and if files exists

  if(opt.kmerWeight != 27)
  {
    throw logic_error("Given kmerWeight value not implemented, select one "
                      "of the following (27), default: 27");
  }

  initialize();

  // TODO: organize following lines, for now first test of workflow


  KmerTranslator *translator = NULL;

  /* get kmer-count and sample files */
  string kmerCountFile = opt.kcFile;
  //vector<string> *sampleList = getFileList(opt.sampleListFile.c_str());
  string seqFile = opt.seqFile;
  //uint8_t readAvgLenList = opt.readAvgLen;

  //for (vector<string>::iterator it = kmerCountList->begin(), jt=readAvgLenList->begin(); it != kmerCountList->end(); ++it, ++jt)
  //{
    /* get dsk kmer-count storage */
    Storage* storage = StorageFactory(STORAGE_HDF5).load(kmerCountFile);
    LOCAL (storage);

    string kmerSizeStr = storage->getGroup("dsk").getProperty ("kmer_size");
    unsigned int kmerSize = atoi(kmerSizeStr.c_str());

    /* just for now, implement other kmerSize later */
    if (kmerSize != 41)
    {
      fprintf(stderr, "kmerSize %u used in hdf5 file %s is not supported yet.\n"
              "For now only dsk output with kmerSize 41 is supported\n",
              kmerSize, kmerCountFile.c_str());
      return EXIT_FAILURE;
    }

    if (translator == NULL)
    {
      translator = new KmerTranslator(kmerSize, opt.kmerWeight);
    }
    //TODO: new translator if (translator->getSpan() != kmerSize)

    unsigned long avgLen = opt.readAvgLen;//stol(*jt);
   
    //TODO: handle mathematical cases 
    float corrFactor = (avgLen==kmerSize)?1:(float)(avgLen)/(avgLen-kmerSize+1);

    fprintf(stderr, "preprocessing...\n");
    /* build lookuptale */
    Lookuptable* lookuptable = buildLookuptable(*storage, *translator, 0, corrFactor);
    if (lookuptable == NULL)
    {
      fprintf(stderr,"Generating lookuptable based on %s failed\n",
              kmerCountFile.c_str());
      return EXIT_FAILURE;
    }
    fprintf(stderr, "finished build lookuptable\n");
    string filename = kmerCountFile;
    size_t lastdot = filename.find_last_of(".");
    if (lastdot != std::string::npos)
      filename=filename.substr(0, lastdot);
    string resultFilename=(string("coverage.")+string(basename(filename.c_str()))+string(".txt"));
    if (opt.threads == 1)
    {
      /*process_sampleList(sampleList, resultFilename, lookuptable, translator,
                         writePopulationCoverage);
      */
      int retval = processSeqFile(seqFile,
                                  resultFilename.c_str(),
                                  *lookuptable,
                                  *translator,
                                  writePopulationCoverage);
      if (retval != EXIT_SUCCESS)
      {
        std::cerr << "ERROR processing read file " << seqFile << std::endl;
      }
    }
    else
    {
      std::cerr << "ERROR: NOT IMPLEMENTED YET" <<std::endl;
      return EXIT_FAILURE;
      /*fprintf(stderr, "start process sampleList\n");
      process_sampleList_threads(sampleList, resultFilename,
                                 lookuptable, translator, opt.threads);*/
    }

    delete lookuptable;
  //}

  //delete sampleList;
  //delete kmerCountList;
  delete translator;

  return EXIT_SUCCESS;
}
