/*
 * Copyright (C) 2018 Annika Seidel <annika.seidel@mpibpc.mpg.de>
 *
 */

#include <cstdio>
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
#include "preprocessing.h"
#include "processReads.h"
#include "filehandling.h"

std::mutex sampleList_mutex;
std::chrono::milliseconds interval(100);

vector<string> *getFileList(const char *fileListFilename)
{
  vector<string> *fileList = new vector<string>();
  string line;
  ifstream fp;
  fp.open(fileListFilename);
  while (getline (fp,line))
    fileList->push_back(line);
  fp.close();

  return fileList;
}

void thread_runner(int id, vector<string> *sampleList, std::string outdir,
                   Lookuptable *lookuptable, KmerTranslator *translator)
{
  std::string sampleFileName;
  string outfilename;
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
      auto path_end = sampleFileName.find_last_of("/");

      if(path_end != std::string::npos)
      {
        outfilename = sampleFileName.substr(path_end+1);
      }
      else
      {
        outfilename = sampleFileName;
      }
      sampleList->pop_back();
      std::cout << "thread " << id << " ready for action, working on "
                << sampleFileName << " outfilename: " << outfilename
                << std::endl << std::flush;
      sampleList_mutex.unlock();
    }
    processReadFile(sampleFileName, outdir + outfilename,
                    *lookuptable, *translator);

  }
  return;
  //FILE * outfile = openFileOrDie(outdir + sampleFileName, "w");
  //fprintf(outfile, "thread %d working on it\n", id);

}

void process_sampleList_threads(vector<string> *sampleList,
                                std::string resultFileName,
                                Lookuptable *lookuptable,
                                KmerTranslator *translator,
                                int num_threads)
{
  vector<string> *cp_sampleList = new vector<string>(*sampleList);
  std::cout << "processing with " << num_threads << " threads" << std::endl
            << std::flush;
  string outdir = "tmp/";
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
  }
}

void process_sampleList(vector<string> *sampleList, std::string resultFileName,
                        Lookuptable *lookuptable,\
                        KmerTranslator *translator)
{
  for (vector<string>::iterator sampleIt = sampleList->begin() ; sampleIt != sampleList->end(); ++sampleIt)
  {

    int retval = processReadFile(sampleIt->c_str(),
                                 resultFileName.c_str(),
                                 *lookuptable,
                                 *translator);
    if (retval != EXIT_SUCCESS)
    {
      std::cerr << "ERROR processing read file " << *sampleIt << std::endl;
    }
  }
}

int pcoverage(int argc, const char **argv, const ToolInfo* tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);
  printf("threads %u\n", opt.threads);
  printf("sampleList: %s\n", opt.sampleListFile.c_str());
  printf("kmerCountList: %s\n", opt.kmerCountListFile.c_str());
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
  vector<string> *kmerCountList = getFileList(opt.kmerCountListFile.c_str());
  vector<string> *sampleList = getFileList(opt.sampleListFile.c_str());
  vector<string> *readAvgLenList = getFileList(opt.readAvgLenFile.c_str());

  for (vector<string>::iterator it = kmerCountList->begin(), jt=readAvgLenList->begin(); it != kmerCountList->end(); ++it, ++jt)
  {
    /* get dsk kmer-count storage */
    Storage* storage = StorageFactory(STORAGE_HDF5).load(*it);
    LOCAL (storage);

    string kmerSizeStr = storage->getGroup("dsk").getProperty ("kmer_size");
    unsigned int kmerSize = atoi(kmerSizeStr.c_str());

    /* just for now, implement other kmerSize later */
    if (kmerSize != 41)
    {
      fprintf(stderr, "kmerSize %u used in hdf5 file %s is not supported yet.\n"
              "For now only dsk output with kmerSize 41 is supported\n",
              kmerSize, it->c_str());
      return EXIT_FAILURE;
    }

    if (translator == NULL)
    {
      translator = new KmerTranslator(kmerSize, opt.kmerWeight);
    }
    //TODO: new translator if (translator->getSpan() != kmerSize)

    unsigned long avgLen = stol(*jt);
    float corrFactor = kmerSize/(avgLen-kmerSize+1);

    /* build lookuptale */
    Lookuptable* lookuptable = buildLookuptable(*storage, *translator, 0, corrFactor);

    std::string resultFileName = "populationCoverages";
    if (opt.threads == 1)
    {
      process_sampleList(sampleList, resultFileName, lookuptable, translator);
    }
    else
    {
      process_sampleList_threads(sampleList, resultFileName,
                                 lookuptable, translator, opt.threads);
    }
    //TODO: correction factor, check retval

    delete lookuptable;
  }

  delete sampleList;
  delete kmerCountList;
  delete translator;

  return EXIT_SUCCESS;
}
