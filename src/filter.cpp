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

#define FREQUENT_COUNT_WINDOWSIZE 10
#define FREQUENT_COUNT_BANDWIDTH 0.2

typedef struct {
  FILE *filterReads;
  FILE *cleanedReads;
  unsigned int dropLevel1;
  std::vector<float> dropLevel2;
  std::vector<uint32_t> shrinkedPlainPositions;
  bool maskOnlyDropEdges;
} FilterArgs;

int filterProcessor(CountProfile &countprofile, void *filterargs)
{

  FilterArgs *currFilterArgs = (FilterArgs *) filterargs;
  SequenceInfo *seqinfo = countprofile.getSeqInfo();

  /* estimate coverage value */
  unsigned int covEst;
  if(currFilterArgs->shrinkedPlainPositions.size() == 0)
    covEst = countprofile.calcMedian();
  else {
    //covEst = countprofile.calcMedian(currFilterArgs->shrinkedPlainPositions);
    covEst = countprofile.calcXquantile(0.67, currFilterArgs->shrinkedPlainPositions);
  }

  Info(Info::DEBUG) << seqinfo->name << "\t" << covEst << "\n";

  uint32_t *maxProfile = countprofile.maximize();

  //TODO: wobbly

  unsigned int mincount = ((FilterArgs *) filterargs)->dropLevel1;
  std::vector<float> percDropLevels = ((FilterArgs *) filterargs)->dropLevel2;

  if(covEst < std::max((unsigned int)10, 2*mincount))
    return 0; //TODO

  bool filter = countprofile.checkForSpuriousTransitionDrops(maxProfile, mincount, currFilterArgs->maskOnlyDropEdges);

  if(!filter) {

    for (float p : percDropLevels) {
      unsigned int dropLevelCriterion = p * covEst;

      if (dropLevelCriterion < 1)
        continue;

      filter = countprofile.checkForSpuriousTransitionDrops(maxProfile, dropLevelCriterion, currFilterArgs->maskOnlyDropEdges);

      if (filter)
        break;
    }
  }

  if (filter)
    sequenceInfo2FileEntry(seqinfo, ((FilterArgs *) filterargs)->filterReads);
  else
    sequenceInfo2FileEntry(seqinfo, ((FilterArgs *) filterargs)->cleanedReads);
}

typedef struct{
  unsigned int numSequences;
  std::vector<uint32_t> summedCounts;
} ProfileStatistic;

int sumCounts(CountProfile &countprofile, void *stat)
{

  ((ProfileStatistic *) stat)->numSequences += 1;
  countprofile.addCountPerPosition(((ProfileStatistic *) stat)->summedCounts);

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
        Info(Info::ERROR) << "ERROR: found invalid argument for " << opt.OP_DROP_LEVEL2.name << " (valid range: 0-0.5)\n";
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

  vector<uint32_t> shrinkedPlainPositions{};
  if (opt.OP_ALIGNED.isSet) {

    Info(Info::INFO) << "try to optimize coverage estimation\n";

    ProfileStatistic stat = {0, std::vector<uint32_t>{}};
    processSeqFile(seqFile, lookuptable, translator, sumCounts, &stat);

    vector<uint32_t> summedCountsSorted = stat.summedCounts;
    sort(summedCountsSorted.begin(), summedCountsSorted.end());

    unsigned int windowSize = FREQUENT_COUNT_WINDOWSIZE;
    double maxDensity = 0, maxDenseCount = 0;
    for (size_t idx = 0; idx < summedCountsSorted.size() - windowSize; idx++) {
      double density =
        (double) windowSize / ((double) (summedCountsSorted[idx + windowSize] - summedCountsSorted[idx]));
      if (density > maxDensity) {
        maxDensity = density;
        maxDenseCount = (summedCountsSorted[idx + windowSize] + summedCountsSorted[idx]) / 2;
      }
    }

    double bandwidth = FREQUENT_COUNT_BANDWIDTH;
    unsigned int sharedCountThr = maxDenseCount + bandwidth * maxDenseCount + 1;

    for (size_t idx = 0; idx < stat.summedCounts.size(); idx++) {
      if (stat.summedCounts[idx] <= sharedCountThr)
        shrinkedPlainPositions.push_back(idx);
    }

    if (shrinkedPlainPositions.size() < 2 * translator->getSpan()) {
      Info(Info::WARNING) << "WARNING: the optimization for coverage estimation (set by --aligned paramter) can not be " \
                           "performed, because less than 2*k positions have counts smaller than the calculated threshold "\
                           "for shared counts. The average profile is probably to spiky to optimize the coverage, " \
                           "estimation. ";
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
    outprefix = getFilename(seqFile);
  string ext = getFileExtension(seqFile);
  FilterArgs filterargs = {openFileOrDie(outprefix + ".spurious" + ext, "w"),
                           openFileOrDie(outprefix + ".cleaned" + ext, "w"),
                           opt.dropLevel1, opt.dropLevel2, shrinkedPlainPositions, maskOnlyDropEdges};

  processSeqFile(seqFile, lookuptable, translator, filterProcessor, &filterargs);

  fclose(filterargs.filterReads);
  fclose(filterargs.cleanedReads);

  delete lookuptable;
  delete translator;

  return EXIT_SUCCESS;
}


