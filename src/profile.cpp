// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>

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
  if (!(seqinfo->comment.empty())) {
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
  opt.printParameterSettings(*tool);

  initialize();
  KmerTranslator *translator = new KmerTranslator(opt.spacedKmerPattern);

  vector<std::string> readFilenames;

  if (opt.OP_READS.isSet)
    readFilenames.push_back(opt.reads);
  if (opt.OP_FORWARD_READS.isSet)
    readFilenames.push_back(opt.forwardReads);
  if (opt.OP_REVERSE_READS.isSet)
    readFilenames.push_back(opt.reverseReads);

  Info(Info::INFO) << "Step 1: Generate lookuptable...\n";
  LookupTableBase *lookuptable;

  // use precomputed counts and fill lookuptable
  if (opt.OP_COUNT_FILE.isSet) {
    string countFile = opt.countFile;
    lookuptable = buildLookuptable(countFile, opt.countMode, *translator, 0);
  } else { // count k-mers itself and fill hash-lookuptable
    lookuptable = buildHashTable(readFilenames, *translator);
  }

  if (lookuptable == NULL) {

    Info(Info::ERROR) << "ERROR: Generating lookuptable failed\n";
    return EXIT_FAILURE;
  }

  ProfileArgs args = {NULL};

  int exit_code = 0;
  Info(Info::INFO) << "Step 2: Print profiles...\n";
  if (opt.threads == 1) {
    if (!opt.reads.empty()) {
      string outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.reads);
      args.profileFile = openFileOrDie(opt.outdir + outprefix + ".profile.reads.txt", "w");
      exit_code = processReads(opt.reads, lookuptable, translator, showProfile, &args, opt.skip);
      fclose(args.profileFile);
    }
    if(exit_code == 0 && !opt.forwardReads.empty()) {
      string outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.forwardReads);
      args.profileFile = openFileOrDie(opt.outdir + outprefix + ".profile.1.txt", "w");
      exit_code = processReads(opt.forwardReads, lookuptable, translator, showProfile, &args, opt.skip);
      fclose(args.profileFile);
    }
    if(exit_code == 0 && !opt.reverseReads.empty()) {
      string outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.reverseReads);
      args.profileFile = openFileOrDie(opt.outdir + outprefix + ".profile.2.txt", "w");
      exit_code = processReads(opt.reverseReads, lookuptable, translator, showProfile, &args, opt.skip);
      fclose(args.profileFile);
    }
  } else {
    // not implemented yet
    assert(false); //should never reach that line
    exit_code = -1;
  }

  opt.deleteInstance();
  delete lookuptable;
  delete translator;

  return exit_code;
}

