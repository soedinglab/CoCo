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

#define FREQUENT_COUNT_WINDOWSIZE 10
#define FREQUENT_COUNT_BANDWIDTH 0.2

typedef struct {
  FILE *filterReads1;
  FILE *filterReads2;
  double threshold;
  /*FILE *cleanedReads;
  FILE *wobblyReads;
  double dropLevel1;
  double dropLevel2;
  std::vector<uint32_t> shrinkedPlainPositions;
  bool maskOnlyDropEdges;*/
} FilterArgs;

int filterProcessor(CountProfile &countprofile, void *filterargs, bool skip)
{
  if(skip) {
    sequenceInfo2FileEntry(countprofile.getSeqInfo(), ((FilterArgs *) filterargs)->filterReads1);
    return 0;
  }

  FilterArgs *currFilterArgs = (FilterArgs *) filterargs;

  bool status = countprofile.checkForSpuriousTransitionDropsWithWindowNew(currFilterArgs->threshold);

  if (status == false)
    sequenceInfo2FileEntry(countprofile.getSeqInfo(), ((FilterArgs *) filterargs)->filterReads1);

  return 0;
}

int filterProcessorPaired(CountProfile &r1_countprofile, CountProfile &r2_countprofile, void *filterargs, bool skip)
{

  if(skip) {
    sequenceInfo2FileEntry(r1_countprofile.getSeqInfo(), ((FilterArgs *) filterargs)->filterReads1);
    sequenceInfo2FileEntry(r2_countprofile.getSeqInfo(), ((FilterArgs *) filterargs)->filterReads2);
    return 0;
  }

  FilterArgs *currFilterArgs = (FilterArgs *) filterargs;

  bool status1, status2;
  status1 = r1_countprofile.checkForSpuriousTransitionDropsWithWindowNew(currFilterArgs->threshold);
  if (status1 == false)
   status2 = r2_countprofile.checkForSpuriousTransitionDropsWithWindowNew(currFilterArgs->threshold);

  if (status1 == false && status2 == false) {
    // no spurious transition drop was found in both paired reads
    sequenceInfo2FileEntry(r1_countprofile.getSeqInfo(), ((FilterArgs *) filterargs)->filterReads1);
    sequenceInfo2FileEntry(r2_countprofile.getSeqInfo(), ((FilterArgs *) filterargs)->filterReads2);
  }
  return 0;
}

typedef struct{
  uint64_t numSequences;
  uint32_t minProfileLen;
  uint32_t maxProfileLen;
  std::vector<uint64_t> summedCounts;
} ProfileStatistic;

int sumCounts(CountProfile &countprofile, void *stat, bool skip)
{
  //TODO: skip
  unsigned int len = countprofile.getProfileLen();
  if (len > ((ProfileStatistic *) stat)->maxProfileLen) {
      ((ProfileStatistic *) stat)->maxProfileLen = len;
  }
  if (len < ((ProfileStatistic *) stat)->minProfileLen ||
      ((ProfileStatistic *) stat)->minProfileLen == 0 ) {
      ((ProfileStatistic *) stat)->minProfileLen = len;
  }

  ((ProfileStatistic *) stat)->numSequences += 1;
  countprofile.addCountPerPosition(((ProfileStatistic *) stat)->summedCounts);

  return 0;
}

int filter(int argc, const char **argv, const Command *tool)
{
  Options &opt = Options::getInstance();
  // overwrite threshold parameter
  opt.OP_THRESHOLD.description="percent drop threshold between two successive windows";
  opt.threshold = 0.1;

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

  SeqInfoMode mode = AUTO;
  for (string readFilename:readFilenames){
    SeqInfoMode currMode = getSeqMode(readFilename);
    if (mode == AUTO)
      mode = currMode;
    else if (currMode != mode) {
      Info(Info::ERROR) << "ERROR: read files have inconsistent file formats\n";
      return EXIT_FAILURE;
    }
  }
  string ext = mode==FASTA?".fa":".fq";

  Info(Info::INFO) << "Step 1: Generate lookuptable...\n";
  LookupTableBase *lookuptable;

  // use precomputed counts and fill lookuptable
  if (opt.OP_COUNT_FILE.isSet) {
    string countFile = opt.countFile;
    lookuptable = buildLookuptable(countFile, opt.countMode, *translator, 0);
    //TODO: change mincount if correction work properly
  } else { // count k-mers itself and fill hash-lookuptable
    lookuptable = buildHashTable(readFilenames, *translator);
  }

  if (lookuptable == NULL) {

    Info(Info::ERROR) << "ERROR: Generating lookuptable failed\n";
    return EXIT_FAILURE;
  }

  FilterArgs args;
  args = {NULL, NULL, opt.threshold};

  int exit_code = 0;
  Info(Info::INFO) << "Step 2: Filter chimeric reads...\n";
  if (!opt.reads.empty()) {
    string outprefix = opt.OP_OUTPREFIX.isSet?opt.outprefix:getFilename(opt.reads);
    args.filterReads1 = openFileOrDie(opt.outdir + outprefix + ".filter.reads" + ext, "w");
    exit_code = processReads(opt.reads, lookuptable, translator, filterProcessor, &args, opt.skip);
    fclose(args.filterReads1);
  }
  if(exit_code == 0 && !opt.forwardReads.empty() && !opt.reverseReads.empty()) {
    string outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.forwardReads);
    args.filterReads1 = openFileOrDie(opt.outdir + outprefix + ".filter.1" + ext, "w");
    outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.reverseReads);
    args.filterReads2 = openFileOrDie(opt.outdir + outprefix + ".filter.2" + ext, "w");
    exit_code = processPairedReads(opt.forwardReads, opt.reverseReads, lookuptable, translator, filterProcessorPaired, &args, opt.skip);
    fclose(args.filterReads1);
    fclose(args.filterReads2);
  }


  Options::deleteInstance();
  delete lookuptable;
  delete translator;

  return exit_code;
}


