// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <cstdio>
#include <stdio.h>
#include <fcntl.h>
#include <stdexcept>
#include <string>

#include <fstream>
#include <iostream>
#include <sys/stat.h>

#include "Command.h"
#include "Info.h"
#include "Options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "preprocessing.h"
#include "runner.h"
#include "filehandling.h"
#include "mask_permuter.h"

#define FREQUENT_COUNT_WINDOWSIZE 10
#define FREQUENT_COUNT_BANDWIDTH 0.2

typedef struct {
  bool dryRun;
  FILE *correctedReadsFasta;
  FILE *correctedReads;
} CorrectorArgs;

int correctionProcessor(CountProfile &countprofile, void *args)
{

  CorrectorArgs *currArgs = (CorrectorArgs *) args;
  SequenceInfo *seqinfo = countprofile.getSeqInfo();

  /* estimate coverage value */
  unsigned int covEst = countprofile.calcMedian();
  Info(Info::DEBUG) << seqinfo->name << "\t" << covEst << "\n";

  /* maximize count profile */
  uint32_t *maxProfile = countprofile.maximize();


  bool changed = false;
  if(covEst > SIGNIFICANT_LEVEL_DIFF)//TODO
      changed = countprofile.correction(maxProfile, covEst, currArgs->dryRun);
  delete[] maxProfile;
  if (currArgs->dryRun) {
    if (changed) {
      fwrite(seqinfo->name.c_str(), sizeof(char), seqinfo->name.size(), currArgs->correctedReads);
      fwrite("\n", sizeof(char), 1, currArgs->correctedReads);
    }
  } else {

    sequenceInfo2FileEntry(seqinfo, currArgs->correctedReadsFasta);
  }
  return 0;
}


int correction(int argc, const char **argv, const Command *tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);


  //check if output file exists
  struct stat buffer;
  string outputfile = opt.outprefix + ".coco_" + tool->cmd + ".txt";
  if (stat(outputfile.c_str(), &buffer) == 0) {
      Info(Info::ERROR) << "ERROR: Outputfile already exists!\n";
      return EXIT_FAILURE;
    }
  //Init mask_permuter class with span 41 and weight 27
  int span = 41;
  int weight = 27;
  if (opt.span != 0){
      span = opt.span;
  }

  if (opt.weight != 0){
      weight = opt.weight;
  }
  mask_permuter mask(span, weight);

  //set permutation range
  unsigned long long int maxPerm = mask.get_permNum();

  //set mask operation range (lexicographically)
  unsigned long long int pmStart = opt.pmstart;
  unsigned long long int pmStop = opt.pmstop;
  unsigned int startPos, stopPos;
  if (pmStop == 0){
      pmStop = maxPerm;
  }
  if (pmStart >= pmStop || pmStart >= maxPerm || pmStop > maxPerm){
      Info(Info::ERROR) << "ERROR: pmstart must be smaller than pmstop"
                           " and should not exceed the maximum of " << maxPerm << "\n";
      return EXIT_FAILURE;
  }

  //if rand is set draw specified number of masks


  opt.dryRun = true; //TODO: change later
  Info(Info::WARNING) << "WARNING: Parameter " << opt.OP_DRY_RUN.display <<
                         " is automatically set to true because the correction step is not implemented yet\n";
  initialize();

  long int perm_count = 0;
  int returnVal;
  unsigned char *msk = new unsigned char[weight];
  std::vector<int> vmsk;
  while(mask.get_next(msk, vmsk)) {
      if((perm_count >= startPos && perm_count <= stopPos) &&
      (perm_count == startPos || perm_count % opt.stepsize == 0)) {
          KmerTranslator *translator = new KmerTranslator();
          //reset translator with span, weight and new mask
          translator->setSW(span, weight);
          translator->setMask(msk);

          //prepare mask for file
          std::string smsk = ">MASK-" + to_string(perm_count) + ":";
          translator->printmask(smsk); //debug check if _mask is set in translator
          for(int i=0; i<weight; i++){
              smsk += std::to_string(vmsk[i]) + " ";
          }

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
          if (opt.OP_OUTPREFIX.isSet) {
              outprefix = opt.outprefix;
          }
          else {
              outprefix = getFilename(seqFile);
          }

          //add mask to file
          FILE * coco_corr = openFileOrDie(outprefix + ".coco_" + tool->cmd + ".txt", "a");
          fwrite(smsk.c_str(), sizeof(char), smsk.length(), coco_corr);
          fwrite("\n", sizeof(char), 1, coco_corr);
          fclose(coco_corr);

          string ext = getFileExtension(seqFile);

          CorrectorArgs args;

          if (opt.dryRun) {
              args = {true,
                      NULL,
                      openFileOrDie(outprefix + ".coco_" + tool->cmd + ".txt", "a")
                      };
              //Info(Info::INFO) << "Perform only a dry run without correction\n";
          } else {
              args = {false,
                      openFileOrDie(outprefix + ".coco_" + tool->cmd + ext, "w"),
                      NULL};
          }

          returnVal = processSeqFile(seqFile, lookuptable, translator, correctionProcessor, &args);

          if (args.correctedReads)
              fclose(args.correctedReads);
          if (args.correctedReadsFasta)
              fclose(args.correctedReadsFasta);

          delete lookuptable;
          //delete translator;
      }
      perm_count++;
  }

  return returnVal;
}

