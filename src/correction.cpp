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
  bool dryRun;
  unsigned int threshold;
  double tolerance;
  CorrectionStatistic *statistic;
  FILE *correctedReadsFasta;
  FILE *errorCandidateReads;
} CorrectorArgs;


int correctionProcessor(CountProfile &countprofile, void *args)
{


  CorrectorArgs *currArgs = (CorrectorArgs *) args;
  SequenceInfo *seqinfo = countprofile.getSeqInfo();

  /* estimate coverage value */
  //unsigned int covEst = countprofile.calcXquantile(0.67);
  //Info(Info::DEBUG) << seqinfo->name << "\t" << covEst << "\n";

  int status = ERROR_FREE;
  CorrectionStatistic *statistic = currArgs->statistic;
  //if(covEst > currArgs->threshold + (unsigned int) (currArgs->tolerance * covEst + 1)) {

    /* maximize count profile */
    uint32_t *maxProfile = countprofile.maximize();
 //sub - indel
    do {
      status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, false, &(statistic->substitution_multikmer), currArgs->dryRun);

      if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
        countprofile.update();
        delete[] maxProfile;
        maxProfile = countprofile.maximize();
        //TODO: update maxProfile instead of generating new one
      }
    } while (!currArgs->dryRun && status == SOME_CORRECTED);

    countprofile.doIndelCorrection(maxProfile, currArgs->threshold, currArgs->tolerance,false, &(statistic->substitution_independent),  &(statistic->insertion),  &(statistic->deletion));
    //TODO: adjust profilelen?
  countprofile.update();
  countprofile.doTrimming(maxProfile, currArgs->threshold, currArgs->tolerance, &(statistic->trimmed));

// sub2 - indel -sub
  /*do {
    status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, true, currArgs->dryRun);

    if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
      countprofile.update();
      delete[] maxProfile;
      maxProfile = countprofile.maximize();
      //TODO: update maxProfile instead of generating new one
    }
  } while (!currArgs->dryRun && status == SOME_CORRECTED);

  countprofile.doIndelCorrection(maxProfile, currArgs->threshold, currArgs->tolerance);
  countprofile.update();
  do {
    status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, false, currArgs->dryRun);

    if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
      countprofile.update();
      delete[] maxProfile;
      maxProfile = countprofile.maximize();
      //TODO: update maxProfile instead of generating new one
    }
  } while (!currArgs->dryRun && status == SOME_CORRECTED);*/

  /*
  // sub2 - indel wo trimming -sub - trimming
  do {
    status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, true, currArgs->dryRun);

    if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
      countprofile.update();
      delete[] maxProfile;
      maxProfile = countprofile.maximize();
      //TODO: update maxProfile instead of generating new one
    }
  } while (!currArgs->dryRun && status == SOME_CORRECTED);

  countprofile.doIndelCorrection(maxProfile, currArgs->threshold, currArgs->tolerance, false);
  countprofile.update();
  do {
    status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, false, currArgs->dryRun);

    if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
      countprofile.update();
      delete[] maxProfile;
      maxProfile = countprofile.maximize();
      //TODO: update maxProfile instead of generating new one
    }
  } while (!currArgs->dryRun && status == SOME_CORRECTED);

  countprofile.doTrimming(maxProfile, currArgs->threshold, currArgs->tolerance);

    //TODO: add trimming strategy for edge errors?
    */
/*
  // sub2 - indel +sub wo trimming -sub - trimming
  do {
    status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, true, &(statistic->substitution_multikmer), currArgs->dryRun);

    if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
      countprofile.update();
      delete[] maxProfile;
      maxProfile = countprofile.maximize();
      //TODO: update maxProfile instead of generating new one
    }
  } while (!currArgs->dryRun && status == SOME_CORRECTED);

  countprofile.doIndelCorrection(maxProfile, currArgs->threshold, currArgs->tolerance, true,  &(statistic->substitution_independent),  &(statistic->insertion),  &(statistic->deletion));
  countprofile.update();
  do {
    status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, false, &(statistic->substitution_singlekmer), currArgs->dryRun);

    if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
      countprofile.update();
      delete[] maxProfile;
      maxProfile = countprofile.maximize();
      //TODO: update maxProfile instead of generating new one
    }
  } while (!currArgs->dryRun && status == SOME_CORRECTED);

  countprofile.doTrimming(maxProfile, currArgs->threshold, currArgs->tolerance, &(statistic->trimmed));*/
/*
  // sub2 - indel +sub wo trimming -sub
  do {
    status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, true, &(statistic->substitution_multikmer), currArgs->dryRun);

    if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
      countprofile.update();
      delete[] maxProfile;
      maxProfile = countprofile.maximize();
      //TODO: update maxProfile instead of generating new one
    }
  } while (!currArgs->dryRun && status == SOME_CORRECTED);

  countprofile.doIndelCorrection(maxProfile, currArgs->threshold, currArgs->tolerance, true, &(statistic->substitution_independent),  &(statistic->insertion),  &(statistic->deletion));
  countprofile.update();
  do {
    status = countprofile.doSubstitutionCorrection(maxProfile, 0, currArgs->threshold, currArgs->tolerance, false, &(statistic->substitution_singlekmer), currArgs->dryRun);

    if (status == SOME_CORRECTED || status == ALL_CORRECTED) {
      countprofile.update();
      delete[] maxProfile;
      maxProfile = countprofile.maximize();
      //TODO: update maxProfile instead of generating new one
    }
  } while (!currArgs->dryRun && status == SOME_CORRECTED);
*/
  // }
  if (currArgs->dryRun) {
    if (status != ERROR_FREE) {
      fwrite(seqinfo->name.c_str(), sizeof(char), seqinfo->name.size(), currArgs->errorCandidateReads);
      fwrite("\n", sizeof(char), 1, currArgs->errorCandidateReads);
    }
  } else {
    sequenceInfo2FileEntry(seqinfo, currArgs->correctedReadsFasta, FASTA);
  }
  return 0;
}


int correction(int argc, const char **argv, const Command *tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

  //TODO: print parameters
  //TODO:check parameter and if files exists


  initialize();
  KmerTranslator *translator = new KmerTranslator(opt.spacedKmerPattern);
  string seqFile = opt.seqFile;


  LookupTableBase *lookuptable;

  // use precomputed counts and fill lookuptable
  if (opt.OP_COUNT_FILE.isSet) {

    string countFile = opt.countFile;
    lookuptable = buildLookuptable(countFile, opt.countMode, *translator, 0);
    //TODO: change mincount if correction work properly
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

  string ext = ".fa"; //getFileExtension(seqFile);
  CorrectorArgs args;
  CorrectionStatistic statistic{0,0,0,0,0,0};

  if (opt.dryRun){
      args = {true, (unsigned int) opt.threshold, opt.tolerance, &statistic, NULL, openFileOrDie(outprefix + ".coco_" + tool->cmd + ".txt", "w")};
      Info(Info::INFO) << "Perform only a dry run without correction\n";
  }
  else {
      args = {false, (unsigned int) opt.threshold, opt.tolerance, &statistic, openFileOrDie(outprefix + ".coco_" + tool->cmd + ext, "w"), NULL};
  }

  FILE *skipReads = openFileOrDie(outprefix + ".coco_" + tool->cmd + "_skipped" + ext, "w");
  int returnVal = processSeqFile(seqFile, lookuptable, translator, correctionProcessor, &args, opt.skip, skipReads );

  // print statistic
  std::cout << "corrections from multiKmer step (substitutions only): " << statistic.substitution_multikmer << std::endl;
  std::cout << "corrections from indel step (total): " << statistic.substitution_independent+statistic.insertion+statistic.deletion << std::endl;
  std::cout << "corrections from indel step (substitutions): " << statistic.substitution_independent << std::endl;
  std::cout << "corrections from indel step (insertions): " << statistic.insertion << std::endl;
  std::cout << "corrections from indel step (deletions): " << statistic.deletion << std::endl;
  std::cout << "corrections from singleKmer step (substitutions only): " << statistic.substitution_singlekmer << std::endl;
  std::cout << "trimmed nucleotides: " << statistic.trimmed << std::endl;


  if (skipReads)
    fclose(skipReads);
  if(args.errorCandidateReads)
    fclose(args.errorCandidateReads);
  if(args.correctedReadsFasta)
    fclose(args.correctedReadsFasta);

  delete lookuptable;
  delete translator;

  return returnVal;
}

