/*
 * Copyright (C) 2018 Annika Seidel <annika.seidel@mpibpc.mpg.de>
 *
 */

#include <cstdio>
#include <fcntl.h>
#include <stdexcept>
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

int showProfile(CountProfile &countprofile, FILE *resultFile)
{
  fprintf(resultFile,"#%s\n", countprofile.getReadName());
  countprofile.showProfile(resultFile);

}

int profile(int argc, const char **argv, const ToolInfo* tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);
  printf("seqFile: %s\n", opt.seqFile.c_str());
  printf("kmerCounFile: %s\n", opt.kcFile.c_str());
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
  //vector<string> *kmerCountList = getFileList(opt.kmerCountListFile.c_str());
  string kmerCountFile = opt.kcFile;
  //vector<string> *sampleList = getFileList(opt.sampleListFile.c_str());
  string seqFile = opt.seqFile;

  //for (vector<string>::iterator it = kmerCountList->begin(); it != kmerCountList->end(); ++it)
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
    //TODO: new translator if (translator->getSpan() != kmerSize

    fprintf(stderr, "preprocessing...\n");
    /* build lookuptale */
    Lookuptable* lookuptable = buildLookuptable(*storage, *translator, 0, 1);
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
    string resultFilename=(string("count_profile.")+string(basename(filename.c_str()))+string(".txt"));

    //process_sampleList(sampleList, resultFilename, lookuptable, translator,showProfile);
    processSeqFile(seqFile, resultFilename, *lookuptable, *translator,showProfile);
    delete lookuptable;
  //}

  //delete sampleList;
  //delete kmerCountList;
  delete translator;

  return EXIT_SUCCESS;
}

