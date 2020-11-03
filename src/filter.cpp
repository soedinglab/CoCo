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
  unsigned int dropLevel1;
  std::vector<float> dropLevel2;
  bool softFilter;
} FilterArgs;

int filterProcessor(CountProfile &countprofile, void *filterargs)
{


  SequenceInfo *seqinfo = countprofile.getSeqInfo();

  //TODO: special median for 16s data
  unsigned int median = countprofile.calcMedian();

  std::cout << seqinfo->name << std::endl;

  uint32_t *maxProfile = countprofile.maximize();

  bool maskOnlyDropEdges = true;
  if (((FilterArgs *) filterargs)->softFilter)
    maskOnlyDropEdges = false;

  //TODO: wobbly

  unsigned int mincount = ((FilterArgs *) filterargs)->dropLevel1;
  std::vector<float> percDropLevels = ((FilterArgs *) filterargs)->dropLevel2;

  if(median < std::max((unsigned int)10, 2*mincount))
    return 0; //TODO

  bool filter = countprofile.checkForSpuriousTransitionDrops(maxProfile, mincount, maskOnlyDropEdges);

  if(!filter) {

    for (float p : percDropLevels) {
      unsigned int dropLevelCriterion = p * median;

      if (dropLevelCriterion < 1)
        continue;

      filter = countprofile.checkForSpuriousTransitionDrops(maxProfile, dropLevelCriterion, maskOnlyDropEdges);

      if (filter)
        break;
    }
  }

  if (filter)
    sequenceInfo2FileEntry(seqinfo, ((FilterArgs *) filterargs)->filterReads);
  else
    sequenceInfo2FileEntry(seqinfo, ((FilterArgs *) filterargs)->cleanedReads);
}

int filter(int argc, const char **argv, const Command *tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

  //TODO: print parameters

  //TODO:check parameter and if files exists
  if (opt.OP_DROP_LEVEL2.isSet) {
    for (float d : opt.dropLevel2){
      if (d >= 0.5){
        Info(Info::ERROR) << "ERROR: found invalid argument for " << opt.OP_DROP_LEVEL2.name << ". (valid range: 0-0.5)\n";
        return EXIT_FAILURE;
      }
    }
  }

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
                           opt.dropLevel1, opt.dropLevel2, opt.softFilter};

  processSeqFile(seqFile, lookuptable, translator, filterProcessor, &filterargs);

  fclose(filterargs.filterReads);
  fclose(filterargs.cleanedReads);

  delete lookuptable;
  delete translator;

  return EXIT_SUCCESS;
}


