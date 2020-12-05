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

  //TODO: print parameters
  //TODO:check parameter and if files exists

  opt.dryRun = true; //TODO: change later
  Info(Info::WARNING) << "WARNING: Parameter " << opt.OP_DRY_RUN.display <<
                         " is automatically set to true because the correction step is not implemented yet\n";

  initialize();
  KmerTranslator *translator = new KmerTranslator();
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

  string ext = getFileExtension(seqFile);

  CorrectorArgs args;

  if (opt.dryRun){
      args = {true, NULL, openFileOrDie(outprefix + ".coco_" + tool->cmd + ".txt", "w")};
      //Info(Info::INFO) << "Perform only a dry run without correction\n";
  }
  else {
      args = {false, openFileOrDie(outprefix + ".coco_" + tool->cmd + ext, "w"), NULL};
  }

  int returnVal = processSeqFile(seqFile, lookuptable, translator, correctionProcessor, &args);

  if(args.correctedReads)
    fclose(args.correctedReads);
  if(args.correctedReadsFasta)
    fclose(args.correctedReadsFasta);

  delete lookuptable;
  delete translator;

  return returnVal;
}
