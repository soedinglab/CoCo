// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>

#include <cstdio>
#include <fcntl.h>
#include <stdexcept>

#include "Command.h"
#include "Info.h"
#include "Options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "preprocessing.h"
#include "runner.h"
#include "filehandling.h"
#include "Lookuptable.h"

int counts2flat(int argc, const char **argv, const Command *tool) {

  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);
  opt.printParameterSettings(*tool);

  initialize();
  KmerTranslator *translator = new KmerTranslator(opt.spacedKmerPattern);
  string reads = opt.reads;

  // use precomputed counts and fill lookuptable
  Info(Info::INFO) << "Step 1: Generate lookuptable...\n";
  LookupTableBase *lookuptable;
  string countFile;

  if (opt.OP_COUNT_FILE.isSet) {
    countFile = opt.countFile;
    lookuptable = buildLookuptable(countFile, opt.countMode, *translator, 0);
  } else {
    Info(Info::ERROR) << "ERROR: Missing count file\n";
    return EXIT_FAILURE;
  }

  if (lookuptable == NULL) {

    Info(Info::ERROR) << "ERROR: Generating lookuptable failed\n";
    return EXIT_FAILURE;
  }

  string outprefix;
  if (opt.OP_OUTPREFIX.isSet)
    outprefix = opt.outprefix;
  else
    outprefix = getFilename(countFile);

  FILE* fp = openFileOrDie(opt.outdir + outprefix + ".counts2flat.tsv", "w");
  lookuptable->iterateOverAll(fp);
  fclose(fp);

  Options::deleteInstance();
  delete lookuptable;
  delete translator;

  return EXIT_SUCCESS;
}