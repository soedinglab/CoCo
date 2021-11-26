// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include <cstdio>
#include <fcntl.h>
#include <stdexcept>
#include <thread>
#include <vector>

#include "Command.h"
#include "Options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "Lookuptable.h"
#include "CountProfile.h"
#include "preprocessing.h"
#include "runner.h"
#include "filehandling.h"
#include "Info.h"

typedef struct {
  FILE *abundanceFile;
} AbundanceEstimatorArgs;


int abundanceEstimatationProcessor(CountProfile &countprofile, void *abundanceargs) {
  // estimate abundance value
  // TODO
  return 0;
}


int abundanceEstimator(int argc, const char **argv, const Command *tool) {

  Info(Info::ERROR) << " Not implemented yet\n";
  return EXIT_SUCCESS;

  int exit_code = 0;
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

  //TODO: print parameters

  // TODO:check parameter and if files exists

  initialize();
  KmerTranslator *translator = new KmerTranslator(opt.spacedKmerPattern);
  string reads = opt.reads;


  LookupTableBase *lookuptable;

  // use precomputed counts and fill lookuptable
  if (opt.OP_COUNT_FILE.isSet) {

    string countFile = opt.countFile;
    lookuptable = buildLookuptable(countFile, opt.countMode, *translator, 0);
  } else { // count k-mers itself and fill hash-lookuptable

    lookuptable = buildHashTable(reads, *translator);
  }

  if (lookuptable == NULL) {

    Info(Info::ERROR) <<"Generating lookuptable failed\n";
    return EXIT_FAILURE;
  }

  AbundanceEstimatorArgs abundanceargs = {openFileOrDie("abundance", "w")};

  if (opt.threads == 1) {
    exit_code = processReads(reads, lookuptable, translator, abundanceEstimatationProcessor, &abundanceargs, opt.skip,
                             NULL);
    if (exit_code != 0) {
      Info(Info::ERROR) << "ERROR processing sequence file " << reads << "\n";
    }
  }



  fclose(abundanceargs.abundanceFile);
  delete lookuptable;
  delete translator;

  return exit_code;
}
