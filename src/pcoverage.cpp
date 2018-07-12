/*
 * Copyright (C) 2018 Annika Seidel <annika.seidel@mpibpc.mpg.de>
 *
 */

#include <cstdio>
#include <stdexcept>
#include <vector>
#include <gatb/gatb_core.hpp>

#include "ToolInfo.h"
#include "options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "Lookuptable.h"
#include "preprocessing.h"
#include "processReads.h"

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

int pcoverage(int argc, const char **argv, const ToolInfo* tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

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
    float corrFactor = (float)(avgLen)/(avgLen-kmerSize+1);

    /* build lookuptale */
    Lookuptable* lookuptable = buildLookuptable(*storage, *translator, 0, corrFactor);
    if(lookuptable != NULL)
    {
      for (vector<string>::iterator sampleIt = sampleList->begin() ; sampleIt != sampleList->end(); ++sampleIt)
      {
        string filename = *it;
        size_t lastdot = filename.find_last_of(".");
        if (lastdot != std::string::npos)
          filename=filename.substr(0, lastdot);
        string resultFilename=(string("coverage.")+string(basename(filename.c_str()))+string(".txt"));
        //std::cout <<test<<std::endl;
        //const char* resultFilename = test.c_str();//(string("coverage.")+*it+string(".txt")).c_str();
        int retval = processReadFile(sampleIt->c_str(),
                                     resultFilename.c_str(),
                                     *lookuptable,
                                     *translator);
      }

      //TODO: check retval
    }
    else
    {
      fprintf(stderr,"Generating lookuptable based on %s failed\n",
              it->c_str());
      return EXIT_FAILURE;
    }

    delete lookuptable;
  }

  delete sampleList;
  delete kmerCountList;
  delete translator;

  return EXIT_SUCCESS;
}
