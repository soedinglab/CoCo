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


typedef struct{
  unsigned int substitution_multikmer;
  unsigned int substitution_singlekmer;
  unsigned int substitution_independent;
  unsigned int insertion;
  unsigned int deletion;
  unsigned int trimmed;
}CorrectionStatistic;

typedef struct {
  unsigned int threshold;
  double tolerance;
  int maxCorrNum;
  int maxTrimLen;
  bool updateLookup;
  CorrectionStatistic *statistic;
  FILE *correctedReads;
} CorrectorArgs;


int correctionProcessor(CountProfile &countprofile, void *args)
{

  CorrectorArgs *currArgs = (CorrectorArgs *) args;
  SequenceInfo *seqinfo = countprofile.getSeqInfo();
  string sequence = seqinfo->seq;

  /* estimate coverage value */
  //unsigned int covEst = countprofile.calcXquantile(0.67);
  //Info(Info::DEBUG) << seqinfo->name << "\t" << covEst << "\n";

  bool updateLookup = currArgs->updateLookup;
  int status = ERROR_FREE;
  CorrectionStatistic statistic = CorrectionStatistic{ 0,0,0,0,0,0}; //currArgs->statistic;
  //if(covEst > currArgs->threshold + (unsigned int) (currArgs->tolerance * covEst + 1)) {

  /* maximize count profile */
  uint32_t *maxProfile = countprofile.maximize();

  /* substitution correction using firstLastUniqueKmerCorrectionStrategy */
  do {
    status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, true, &(statistic.substitution_multikmer));

    if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
      countprofile.update(updateLookup);
      delete[] maxProfile;
      maxProfile = countprofile.maximize();
    }
  } while (status == SOME_CORRECTED);

  /* indel correction and edge substitution correction */
  bool changed = countprofile.doIndelCorrection(maxProfile, currArgs->threshold, currArgs->tolerance, true, &(statistic.substitution_independent),  &(statistic.insertion),  &(statistic.deletion));
  if(changed) {
    countprofile.update(updateLookup);
    delete[] maxProfile;
    maxProfile = countprofile.maximize();
  }

  /* second round of substitution correction without forcing the use of multiple kmers*/
  if (status != ERROR_FREE) {
    do {
      status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, false,
                                                     &(statistic.substitution_singlekmer));

      if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
        countprofile.update(updateLookup);
        delete[] maxProfile;
        maxProfile = countprofile.maximize();
      }
    } while (status == SOME_CORRECTED);
  }

  /* trimming */
  if(currArgs->maxTrimLen > 0) {
    if(countprofile.doTrimming(maxProfile, currArgs->threshold, currArgs->tolerance, currArgs->maxTrimLen,
                            &(statistic.trimmed)))
      countprofile.update(updateLookup);
  }

  if(currArgs->maxCorrNum > 0 && (int)(statistic.substitution_multikmer+ statistic.substitution_singlekmer + statistic.substitution_independent +
     statistic.insertion + statistic.deletion + statistic.trimmed) > currArgs->maxCorrNum){
    //revert corrections because too many corrections were applied (strain flipping?)

    seqinfo->seq=sequence;
    //countprofile.update(updateLookup);
  } else{

    // update global count statistic
    currArgs->statistic->substitution_multikmer += statistic.substitution_multikmer;
    currArgs->statistic->substitution_singlekmer += statistic.substitution_singlekmer;
    currArgs->statistic->substitution_independent += statistic.substitution_independent;
    currArgs->statistic->insertion += statistic.insertion;
    currArgs->statistic->deletion += statistic.deletion;
    currArgs->statistic->trimmed += statistic.trimmed;
  }

  sequenceInfo2FileEntry(seqinfo, currArgs->correctedReads, AUTO);


  delete[] maxProfile;
  return 0;
}


int correction(int argc, const char **argv, const Command *tool)
{
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

  CorrectorArgs args;
  CorrectionStatistic statistic{0,0,0,0,0,0};
  args = {(unsigned int) opt.threshold, opt.tolerance, opt.maxCorrNum, opt.maxTrimLen, opt.updateLookup, &statistic, NULL};

  FILE *skipReads = openFileOrDie((opt.OP_OUTPREFIX.isSet?opt.outprefix:std::string("")) + ".skipped" + ext, "w");
  int returnVal=0;
  Info(Info::INFO) << "Step 2: sequencing error correction...\n";
  if (!opt.reads.empty()) {
    string outprefix = opt.OP_OUTPREFIX.isSet?opt.outprefix:getFilename(opt.reads);
    //args.correctedReads = openFileOrDie(outprefix + ".coco_" + tool->cmd + ".reads" + ext, "w");
    args.correctedReads = openFileOrDie(outprefix + ".corr.reads" + ext, "w");
    returnVal = processReads(opt.reads, lookuptable, translator, correctionProcessor, &args, opt.skip, skipReads);
    fclose(args.correctedReads);
  }
  if(!opt.forwardReads.empty()) {
    string outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.forwardReads);
    args.correctedReads = openFileOrDie(outprefix + ".corr.1" + ext, "w");
    returnVal = processReads(opt.forwardReads, lookuptable, translator, correctionProcessor, &args, opt.skip,
                             skipReads);
    fclose(args.correctedReads);
  }
  if(!opt.reverseReads.empty()) {
    string outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.reverseReads);
    args.correctedReads = openFileOrDie(outprefix + ".corr.2" + ext, "w");
    returnVal = processReads(opt.reverseReads, lookuptable, translator, correctionProcessor, &args, opt.skip,
                             skipReads);
    fclose(args.correctedReads);
  }

  if (skipReads)
    fclose(skipReads);

  // print statistic
  std::cout << "### COCO ERROR CORRECTION STATISTIC ###" << std::endl;
  std::cout << "corrections from multiKmer step (substitutions only): " << statistic.substitution_multikmer << std::endl;
  std::cout << "corrections from indel step (total): " << statistic.substitution_independent+statistic.insertion+statistic.deletion << std::endl;
  std::cout << "corrections from indel step (substitutions): " << statistic.substitution_independent << std::endl;
  std::cout << "corrections from indel step (insertions): " << statistic.insertion << std::endl;
  std::cout << "corrections from indel step (deletions): " << statistic.deletion << std::endl;
  std::cout << "corrections from singleKmer step (substitutions only): " << statistic.substitution_singlekmer << std::endl;
  std::cout << "trimmed nucleotides: " << statistic.trimmed << std::endl;

  delete lookuptable;
  delete translator;

  return returnVal;
}

