/*
 * Copyright (C) 2018 Annika Seidel <annika.seidel@mpibpc.mpg.de>
 *
 */

#include <cstdio>
#include <stdexcept>
#include <vector>
#include "ToolInfo.h"
#include "options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "Lookuptable.h"
#include "preprocessing.h"

int pcoverage(int argc, const char **argv, const ToolInfo* tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

  printf("argument: sampleList: %s\n", opt.sampleListFile.c_str());
  printf("argument: kmerCountList: %s\n", opt.kmerCountListFile.c_str());
  printf("argument: kmerWeight: %u\n", opt.kmerWeight);

  // TODO:check parameter and if files exists

  if(opt.kmerWeight != 27)
  {
    throw logic_error("Given kmerWeight value not implemented, select one "
                      "of the following (27), default: 27");
  }

  initialize();

  // TODO: organize following lines, for now first test of workflow


  KmerTranslator *translator = NULL;

  /* get kmer-count files */
  const char * kmerCountListFile = opt.kmerCountListFile.c_str();
  vector<string> kmerCountList;
  string line;
  ifstream kmerCountListfp;
  kmerCountListfp.open(kmerCountListFile);
  while (getline (kmerCountListfp,line))
    kmerCountList.push_back(line);
  kmerCountListfp.close();


  for (std::vector<string>::iterator it = kmerCountList.begin() ; it != kmerCountList.end(); ++it)
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

    /* build lookuptale */
    Lookuptable* lookuptable = buildLookuptable(*storage, *translator, 0);

    //TODO: algorithm

    delete lookuptable;
  }

  delete translator;

  return EXIT_SUCCESS;
}
