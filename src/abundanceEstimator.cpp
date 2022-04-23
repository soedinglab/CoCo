// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>
#include <cstdio>
#include <fcntl.h>
#include <stdexcept>
#include <vector>

#include "Command.h"
#include "Options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "Lookuptable.h"
#include "CountProfile.h"
#include "preprocessing.h"
#include "runner.h"
#include "filehandling.h"
#include "Info.h"

typedef struct {
  FILE *abundanceFile;
} AbundanceEstimatorArgs;


int abundanceEstimatationProcessor(CountProfile &countprofile, void *args, bool skip) {
  // estimate abundance value as 67% quantile
  double quantile = 0.67;

  AbundanceEstimatorArgs *currArgs = (AbundanceEstimatorArgs *) args;
  SequenceInfo *seqinfo = countprofile.getSeqInfo();

  if (skip) {
    /* sequence is too short for abundance estimation just write sequence to the outputfile without doing anything */
    Info(Info::CDEBUG) << "WARNING: sequence " << seqinfo->name << " is too short, it'll be skipped\n";
    fprintf(currArgs->abundanceFile, "%s\t-\n", (seqinfo->name).c_str());
    return 0;
  }

  unsigned int est = countprofile.calcXquantile(quantile);
  fprintf(currArgs->abundanceFile, "%s\t%u\n", (seqinfo->name).c_str(), est);
  return 0;
}


int abundanceEstimator(int argc, const char **argv, const Command *tool) {

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

  Info(Info::INFO) << "Step 1: Generate lookuptable...\n";
  LookupTableBase *lookuptable;

  // use precomputed counts and fill lookuptable
  if (opt.OP_COUNT_FILE.isSet) {
    string countFile = opt.countFile;
    lookuptable = buildLookuptable(countFile, opt.countMode, *translator, 0);
  } else { // count k-mers itself and fill hash-lookuptable
    lookuptable = buildHashTable(readFilenames, *translator);
  }

  if (lookuptable == NULL) {

    Info(Info::ERROR) << "ERROR: Generating lookuptable failed\n";
    return EXIT_FAILURE;
  }


  AbundanceEstimatorArgs args = {NULL};

  int exit_code = 0;
  Info(Info::INFO) << "Step 2: Abundance estimation...\n";
  if (opt.threads == 1) {
    if (!opt.reads.empty()) {
      string outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.reads);
      args.abundanceFile = openFileOrDie(opt.outdir + outprefix + ".abundance.reads.tsv", "w");
      exit_code = processReads(opt.reads, lookuptable, translator, abundanceEstimatationProcessor, &args, opt.skip);
      fclose(args.abundanceFile);
    }
    if(exit_code == 0 && !opt.forwardReads.empty()) {
      string outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.forwardReads);
      args.abundanceFile = openFileOrDie(opt.outdir + outprefix + ".abundance.1.tsv", "w");
      exit_code = processReads(opt.forwardReads, lookuptable, translator, abundanceEstimatationProcessor, &args, opt.skip);
      fclose(args.abundanceFile);
    }
    if(exit_code == 0 && !opt.reverseReads.empty()) {
      string outprefix = opt.OP_OUTPREFIX.isSet ? opt.outprefix : getFilename(opt.reverseReads);
      args.abundanceFile = openFileOrDie(opt.outdir + outprefix + ".abundance.2.tsv", "w");
      exit_code = processReads(opt.reverseReads, lookuptable, translator, abundanceEstimatationProcessor, &args, opt.skip);
      fclose(args.abundanceFile);
    }
  } else {
    // not implemented yet
    assert(false); //should never reach that line
    exit_code = -1;
  }

  opt.deleteInstance();
  delete lookuptable;
  delete translator;

  return exit_code;
}
