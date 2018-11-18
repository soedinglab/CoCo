/*
 * Copyright (C) 2018 Annika Seidel <annika.seidel@mpibpc.mpg.de>
 *
 */

#include <cstdio>
#include <fcntl.h>
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

std::chrono::milliseconds interval(100);


int writeAbundanceEstimation(CountProfile &countprofile, FILE *resultFile)
{
  /* estimate population coverage */
  float abundanceEstimation = countprofile.calc67quantile() * countprofile.getCorrFactor();
  fprintf(resultFile, "%s\t%f\n", countprofile.getSeqName(), abundanceEstimation );
}


int abundanceEstimator(int argc, const char **argv, const ToolInfo* tool)
{
  int retval = 0;
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);
  printf("threads %u\n", opt.threads);
  printf("seqFile: %s\n", opt.seqFile.c_str());
  printf("kcFile: %s\n", opt.kcFile.c_str());
  printf("kmerWeight: %u\n", opt.kmerWeight);
  
  // TODO:check parameter and if files exists
  // TODO: check resultfile already exists and tmp folder

  if(opt.kmerWeight != 27)
  {
    throw logic_error("Given kmerWeight value not implemented, select one "
                      "of the following (27), default: 27");
  }

  initialize();
  KmerTranslator *translator = NULL;
  Lookuptable* lookuptable = NULL;
  // TODO: organize following lines (preprocessing funct), for now first test of workflow




  /* get kmer-count and sample files */
  string kmerCountFile = opt.kcFile;
  string seqFile = opt.seqFile;


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

  translator = new KmerTranslator(kmerSize, opt.kmerWeight);
  //TODO: check return

  float corrFactor = 1;
  if (opt.readAvgLen != 0)
  {
    unsigned long avgLen = opt.readAvgLen;//stol(*jt); //TODO float
    //TODO: handle mathematical cases
    corrFactor = (avgLen==kmerSize)?1:(float)(avgLen)/(avgLen-kmerSize+1);
  }

  fprintf(stderr, "preprocessing...\n");
  /* build lookuptale */
  lookuptable = buildLookuptable(*storage, *translator, 0, corrFactor);
  if (lookuptable == NULL)
  {
    fprintf(stderr,"Generating lookuptable based on %s failed\n",
            kmerCountFile.c_str());
    return EXIT_FAILURE;
  }
  fprintf(stderr, "finished build lookuptable\n");

  /* process seq file */
  if (retval == 0)
  {

    string resultFile = get_filename(seqFile) + string(".") +
                        get_filename(kmerCountFile) + string(".abundance");

    printf("resultFilename: %s\n", resultFile.c_str());
    if (opt.threads == 1)
    {
      retval = processSeqFile(seqFile,
                              resultFile,
                              lookuptable,
                              translator,
                              writeAbundanceEstimation);
      if (retval != EXIT_SUCCESS)
      {
        std::cerr << "ERROR processing sequence file " << seqFile << std::endl;
      }
    }
    else
    {
      retval = processSeqFileParallel(seqFile,
                                      resultFile,
                                      lookuptable,
                                      translator,
                                      writeAbundanceEstimation,
                                      opt.threads);
    }
  }

  delete lookuptable;
  delete translator;

  return retval;
}
