// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <cstdio>
#include <fcntl.h>
#include <stdexcept>
#include <vector>

#include "Command.h"
#include "options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "CountProfile.h"
#include "preprocessing.h"
#include "runner.h"
#include "filehandling.h"


typedef struct {
  FILE *profileFile;
} ProfileArgs;


int showProfile(CountProfile &countprofile, void *profileargs) {
  fprintf(((ProfileArgs *) profileargs)->profileFile, "#%s\n", countprofile.getSeqName());
  countprofile.showProfile(((ProfileArgs *) profileargs)->profileFile);

}

int profile(int argc, const char **argv, const Command *tool) {

  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

  //TODO: print parameters

  // TODO:check parameter and if files exists

  initialize();
  KmerTranslator *translator = new KmerTranslator();
  string seqFile = opt.seqFile;


  LookupTableBase *lookuptable;

  // use precomputed counts and fill lookuptable
  if (opt.OP_COUNT_FILE.isSet) {

    string countFile = opt.countFile;
    lookuptable = buildLookuptable(countFile, *translator, 0, 1);
  } else { // count k-mers itself and fill hash-lookuptable

    lookuptable = buildHashTable(seqFile, *translator);
  }

  if (lookuptable == NULL) {

    fprintf(stderr, "Generating lookuptablefailed\n");
    return EXIT_FAILURE;
  }

  string outprefix;
  if (opt.OP_OUTPREFIX.isSet)
    outprefix = opt.outprefix;
  else
    outprefix = getFilename(seqFile);

  ProfileArgs profileargs = {openFileOrDie(outprefix + ".profiles", "w")};
  processSeqFile(seqFile, lookuptable, translator, showProfile, &profileargs);
  fclose(profileargs.profileFile);

  return EXIT_SUCCESS;
}

