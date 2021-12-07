// Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>

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
  FILE *filterReads;
  FILE *cleanedReads;
  FILE *wobblyReads;
  double dropLevel1;
  double dropLevel2;
  std::vector<uint32_t> shrinkedPlainPositions;
  bool maskOnlyDropEdges;
} FilterArgs;

int filterProcessor(CountProfile &countprofile, void *filterargs, bool skip)
{
  // TODO: skip logic

  FilterArgs *currFilterArgs = (FilterArgs *) filterargs;
  SequenceInfo *seqinfo = countprofile.getSeqInfo();

  /* estimate coverage value */
  unsigned int covEst;
  if(currFilterArgs->shrinkedPlainPositions.empty()){

    covEst = countprofile.calcMedian();
  } else {

    covEst = countprofile.calcXquantile(0.67, currFilterArgs->shrinkedPlainPositions);
  }

  Info(Info::CDEBUG) << seqinfo->name << "\t" << covEst << "\n";

  uint32_t *maxProfile = countprofile.maximize();



  if(covEst <= SIGNIFICANT_LEVEL_DIFF) {
    delete[] maxProfile;
    sequenceInfo2FileEntry(seqinfo, ((FilterArgs *) filterargs)->wobblyReads);
    return 0;
  }

  float dropLevel1 = ((FilterArgs *) filterargs)->dropLevel1;
  float dropLevel2 = ((FilterArgs *) filterargs)->dropLevel2;
  bool maskOnlyDropEdges = ((FilterArgs *) filterargs)->maskOnlyDropEdges;
  bool filter = countprofile.checkForSpuriousTransitionDropsWithWindow(maxProfile, covEst, dropLevel1, dropLevel2, maskOnlyDropEdges);
  delete[] maxProfile;


  if (filter)
    sequenceInfo2FileEntry(seqinfo, ((FilterArgs *) filterargs)->filterReads);
  else
    sequenceInfo2FileEntry(seqinfo, ((FilterArgs *) filterargs)->cleanedReads);

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
  opt.parseOptions(argc, argv, *tool);

  //TODO: print parameters

  //TODO:check parameter and if files exists
  if (opt.OP_DROP_LEVEL1.isSet && (opt.dropLevel1 < 0 || opt.dropLevel1 > 0.33)){
      Info(Info::ERROR) << "ERROR: found invalid argument for " << opt.OP_DROP_LEVEL1.name << " (valid range: 0-0.33)\n";
      return EXIT_FAILURE;
  }
  if (opt.OP_DROP_LEVEL2.isSet && (opt.dropLevel2 < 0 || opt.dropLevel2 > 0.33)){
      Info(Info::ERROR) << "ERROR: found invalid argument for " << opt.OP_DROP_LEVEL2.name << " (valid range: 0-0.33)\n";
      return EXIT_FAILURE;
  }

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

  vector<uint32_t> shrinkedPlainPositions{};
  if (opt.OP_ALIGNED.isSet) {

    Info(Info::INFO) << "prepare for optimized coverage estimation\n";

    ProfileStatistic stat = {0, 0, 0, std::vector<uint64_t>{}};

    processReads(reads, lookuptable, translator, sumCounts, &stat, true);

    vector<uint64_t> summedCountsSorted = stat.summedCounts;
    summedCountsSorted.resize(stat.minProfileLen);
    sort(summedCountsSorted.begin(), summedCountsSorted.end());

    unsigned int windowSize = FREQUENT_COUNT_WINDOWSIZE;
    double maxDensity = 0, maxDenseCount = 0;
    for (size_t idx = 0; idx < summedCountsSorted.size() - windowSize; idx++) {
      double density =
        (double) windowSize / ((double) (summedCountsSorted[idx + windowSize] - summedCountsSorted[idx] + 1));
      if (density > maxDensity) {
        maxDensity = density;
        maxDenseCount = (summedCountsSorted[idx + windowSize] + summedCountsSorted[idx]) / 2;
      }
    }

    double bandwidth = FREQUENT_COUNT_BANDWIDTH;
    uint64_t sharedCountThr = maxDenseCount + bandwidth * maxDenseCount + 1;

    for (size_t idx = 0; idx < stat.summedCounts.size(); idx++) {
      if (stat.summedCounts[idx] <= sharedCountThr)
        shrinkedPlainPositions.push_back(idx);
    }

    if (shrinkedPlainPositions.size() < 2 * translator->getSpan()) {
      Info(Info::WARNING) << "WARNING: the optimization for coverage estimation (set by --aligned paramter) can not be " \
                           "performed, because less than 2*k positions have counts smaller than the calculated pseudocount "\
                           "for shared counts. The average profile is probably to spiky to optimize the coverage, " \
                           "estimation.\n";
      shrinkedPlainPositions.clear();
      /*shrinkedPlainPositions.reserve(stat.summedCounts.size());
      for (size_t idx = 0; idx < stat.summedCounts.size(); idx++) {
        shrinkedPlainPositions.push_back(idx);
      }*/
    }
  }

  bool maskOnlyDropEdges = true;
  if (opt.softFilter)
    maskOnlyDropEdges = false;

  string outprefix;
  if (opt.OP_OUTPREFIX.isSet)
    outprefix = opt.outprefix;
  else
    outprefix = getFilename(reads);

  string ext = getFileExtension(reads);
  FilterArgs filterargs = {openFileOrDie(outprefix + ".coco_" + tool->cmd + ".spurious" + ext, "w"),
                           openFileOrDie(outprefix + ".coco_" + tool->cmd + ".cleaned" + ext, "w"),
                           openFileOrDie(outprefix + ".coco_" + tool->cmd + ".wobbly" + ext, "w"),
                           opt.dropLevel1, opt.dropLevel2, shrinkedPlainPositions, maskOnlyDropEdges};

  FILE *skipReads = openFileOrDie(outprefix + ".coco_" + tool->cmd + "_skipped" + ext, "w");
  int returnVal = processReads(reads, lookuptable, translator, filterProcessor, &filterargs, opt.skip);

  fclose(skipReads);
  fclose(filterargs.filterReads);
  fclose(filterargs.wobblyReads);
  fclose(filterargs.cleanedReads);

  opt.deleteInstance();
  delete lookuptable;
  delete translator;

  return returnVal;
}


