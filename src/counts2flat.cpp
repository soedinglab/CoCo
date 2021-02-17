// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <cstdio>
#include <stdio.h>
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

  //TODO: print parameters
  //TODO:check parameter and if files exists

  initialize();
  KmerTranslator *translator = new KmerTranslator();
  string seqFile = opt.seqFile;

  // use precomputed counts and fill lookuptable
  string countFile = opt.countFile;
  Lookuptable *lookuptable = (Lookuptable*) buildLookuptable(countFile, opt.countMode, *translator, 0);
  //TODO: change mincount if correction work properly

  if (lookuptable == NULL) {

    Info(Info::ERROR) << "ERROR: Generating lookuptable failed\n";
    return EXIT_FAILURE;
  }

  string outprefix;
  if (opt.OP_OUTPREFIX.isSet)
    outprefix = opt.outprefix;
  else
    outprefix = getFilename(countFile);

  FILE* fp = openFileOrDie(outprefix + ".coco_" + tool->cmd + ".txt", "w");
  lookuptable->iterateOverAll(fp);

  fclose(fp);
  delete lookuptable;
  delete translator;

  return EXIT_SUCCESS;
}