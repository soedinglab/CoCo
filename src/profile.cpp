// Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>

#include <cstdio>
#include <fcntl.h>
#include <stdexcept>
#include <vector>

#include "Command.h"
#include "Info.h"
#include "Options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "CountProfile.h"
#include "preprocessing.h"
#include "runner.h"
#include "filehandling.h"


typedef struct {
  FILE *profileFile;
} ProfileArgs;


int showProfile(CountProfile &countprofile, void *profileargs, bool skip) {

  if (skip)
    return 0;
  FILE *fp = ((ProfileArgs *) profileargs)->profileFile;
  SequenceInfo *seqinfo = countprofile.getSeqInfo();

  fwrite("#", sizeof(char), 1, fp);
  fwrite(seqinfo->name.c_str(), sizeof(char), seqinfo->name.size(), fp);
  if (seqinfo->comment.size() > 0) {
    fwrite(" ", sizeof(char), 1, fp);
    fwrite(seqinfo->comment.c_str(), sizeof(char), seqinfo->comment.size(), fp);
  }
  fwrite("\n", sizeof(char), 1, fp);
  countprofile.showProfile(fp);
  //countprofile.showMaximzedProfile(fp);
  return 0;
}

int profile(int argc, const char **argv, const Command *tool) {

  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

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

    Info(Info::ERROR) << "ERROR: Generating lookuptable failed\n";
    return EXIT_FAILURE;
  }

  string outprefix;
  if (opt.OP_OUTPREFIX.isSet)
    outprefix = opt.outprefix;
  else
    outprefix = getFilename(reads);

  ProfileArgs profileargs = {openFileOrDie(outprefix + ".profiles", "w")};
  processReads(reads, lookuptable, translator, showProfile, &profileargs, opt.skip, NULL);
  fclose(profileargs.profileFile);
  opt.deleteInstance();

  return EXIT_SUCCESS;
}

