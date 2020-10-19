// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <cstdio>
#include <stdio.h>
#include <fcntl.h>
#include <stdexcept>

#include "Command.h"
#include "Info.h"
#include "options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "preprocessing.h"
#include "runner.h"
#include "filehandling.h"

typedef struct {
  FILE *filterReads;
  FILE *cleanedReads;
  int minCount;
} FilterArgs;

int filterProcessor(CountProfile &countprofile, void *filterargs)
{

  SequenceInfo *seqinfo = countprofile.getSeqInfo();
  bool toFilter = countprofile.checkForTransitionDrops(((FilterArgs *) filterargs)->minCount);

  if (toFilter)
    sequenceInfo2FileEntry(seqinfo, ((FilterArgs *) filterargs)->filterReads);
  else
    sequenceInfo2FileEntry(seqinfo, ((FilterArgs *) filterargs)->cleanedReads);
}

int filter(int argc, const char **argv, const Command *tool)
{
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

    Info(Info::ERROR) << "ERROR: Generating lookuptable failed\n";
    return EXIT_FAILURE;
  }

  string outprefix;
  if (opt.OP_OUTPREFIX.isSet)
    outprefix = opt.outprefix;
  else
    outprefix = getFilename(seqFile);
  string ext = getFileExtension(seqFile);
  FilterArgs filterargs = {openFileOrDie(outprefix + ".spurious" + ext, "w"),
                           openFileOrDie(outprefix + ".cleaned" + ext, "w"),
                           opt.minCount};

  processSeqFile(seqFile, lookuptable, translator, filterProcessor, &filterargs);

  fclose(filterargs.filterReads);
  fclose(filterargs.cleanedReads);

  delete lookuptable;
  delete translator;

  return EXIT_SUCCESS;
}


