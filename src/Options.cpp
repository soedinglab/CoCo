// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

#include "Options.h"
#include "Info.h"
#include "util.h"


Options *Options::instance = NULL;

Options::Options() :
  OP_SEQ_FILE(OP_SEQ_FILE_ID, "seqfile", "--seqfile",
               "sequence file (reads or contigs in fasta/fastq format)",
               typeid(std::string), (void *) &seqFile, PROFILE | FILTER | ABUNDANCE_ESTIMATOR | CORRECTOR | CONSENSUS),
  OP_COUNT_FILE(OP_COUNT_FILE_ID, "counts", "--counts",
               "pre computed kmer count file in hdf5 format (dsk output), Note: only supports 41-mers yet",
               typeid(std::string), (void *) &countFile, COUNTS2FLAT),
  OP_OUTPREFIX(OP_OUTPREFIX_ID, "outprefix", "--outprefix",
               "prefix to use for resultfile(s)",
               typeid(std::string), (void *) &outprefix, 0),
  OP_SPACED_KMER_PATTERN(OP_SPACED_KMER_PATTERN_ID,"spaced-pattern",
               "--spaced-pattern", "User-specified spaced k-mer pattern (span must be <=64, 12 <= weight <=32 and symmetric),\n "\
               "default: 11110111111011011101010111011011111101111",
               typeid(std::string), (void*) &spacedKmerPattern, 0),
  OP_SKIP(OP_SKIP_ID,"skip", "--skip", "skip sequences with less than this many k-mers",
               typeid(int), (void*) &skip, 0),
  OP_THRESHOLD(OP_THRESHOLD_ID, "threshold", "--threshold",
               "untrusted count threshold (default: 3)",
               typeid(int), (void *) &threshold, 0),
  OP_TOLERANCE(OP_TOLERANCE_ID, "tolerance", "--tolerance",
               "relative neighborhood count added to threshold value (default: 0.01)",
               typeid(double), (void *) &tolerance, 0),
  OP_MAX_CORR_NUM(OP_MAX_CORR_NUM_ID, "max-corr-num", "--max-corr-num",
               "maximal number of corrections per read, changes are discarded otherwise (default: 0=off)",
               typeid(int), (void *) &maxCorrNum, 0),
  OP_MAX_TRIM_LEN(OP_MAX_TRIM_LEN_ID, "max-trim-len", "--max-trim-len",
               "trim up to this many nucleotides from the beginning/end of reads if no correction was possible (default: 0)",
               typeid(int), (void *) &maxTrimLen, 0),
  OP_UPDATE_LOOKUPTABLE(OP_UPDATE_LOOKUPTABLE_ID, "update-lookup", "--update-lookup",
               "update counts in lookuptable after sequence is corrected (slow down, might help in low coverage regions)",
                typeid(bool), (void *) &updateLookup, 0),
  OP_DROP_LEVEL1(OP_DROP_LEVEL1_ID, "drop-level1", "--drop-level1",
               "local drop criterion (range 0-0.33)",
               typeid(double), (void *) &dropLevel1, 0),
  OP_DROP_LEVEL2(OP_DROP_LEVEL2_ID, "drop-level2", "--drop-level2",
               "global drop criterion (range 0-0.33)",
               typeid(double), (void *) &dropLevel2, 0),
  OP_SOFT_FILTER(OP_SOFT_FILTER_ID, "soft", "--soft",
               "less strict filtering mode due to more strict masking strategy ",
               typeid(bool), (void *) &softFilter, 0),
  OP_ALIGNED(OP_ALIGNED_ID, "aligned", "--aligned",
               "optimize abundance estimation for reads that span the same region (amplicon sequence data)",
               typeid(bool), (void *) &aligned, 0),
  OP_DRY_RUN(OP_DRY_RUN_ID, "dry-run", "--dry-run",
             "perform a trial run that doesn't make any changes", typeid(bool), (void *) &dryRun, 0),
  OP_THREADS(OP_THREADS_ID, "threads", "--threads", "number of threads, not supported yet (default: 1)", typeid(int), (void *) &threads, 0),
  OP_VERBOSE(OP_VERBOSE_ID, "verbose", "--verbose", "verbosity level, 0: quiet 1: Errors, 2: +Warnings, 3: +Info, 4: +Debug, "\
                            "default: 3", typeid(int), (void *) &verbose, 0),
  // expert options
  OP_COUNT_MODE(OP_COUNT_MODE_ID, "count-mode", "--count-mode",
                "way to store counts for concurrent spaced kmers (expert option)\n 0: sum (default)\n 1: maximize",
                typeid(int), (void *) &countMode, 0)

  {
  if (instance) {
    std::cerr << "Parameter instance already exists!\n";
    abort();
  }
  instance = this;

  setDefaults();

  //TODO: threads

  //corrector
  correctionWorkflow.push_back(&OP_SEQ_FILE);
  correctionWorkflow.push_back(&OP_COUNT_FILE);
  correctionWorkflow.push_back(&OP_COUNT_MODE);
  correctionWorkflow.push_back(&OP_OUTPREFIX);
  correctionWorkflow.push_back(&OP_SPACED_KMER_PATTERN);
  correctionWorkflow.push_back(&OP_SKIP);
  correctionWorkflow.push_back(&OP_THRESHOLD);
  correctionWorkflow.push_back(&OP_TOLERANCE);
  correctionWorkflow.push_back(&OP_MAX_CORR_NUM);
  correctionWorkflow.push_back(&OP_MAX_TRIM_LEN);
  correctionWorkflow.push_back(&OP_UPDATE_LOOKUPTABLE);
  correctionWorkflow.push_back(&OP_DRY_RUN);
  correctionWorkflow.push_back(&OP_VERBOSE);

  // filter
  filterWorkflow.push_back(&OP_SEQ_FILE);
  filterWorkflow.push_back(&OP_COUNT_FILE);
  filterWorkflow.push_back(&OP_OUTPREFIX);
  filterWorkflow.push_back(&OP_SPACED_KMER_PATTERN);
  filterWorkflow.push_back(&OP_SKIP);
  filterWorkflow.push_back(&OP_DROP_LEVEL1);
  filterWorkflow.push_back(&OP_DROP_LEVEL2);
  filterWorkflow.push_back(&OP_ALIGNED);
  filterWorkflow.push_back(&OP_SOFT_FILTER);
  filterWorkflow.push_back(&OP_COUNT_MODE);
  filterWorkflow.push_back(&OP_VERBOSE);

  //abundanceEstimator
  abundanceEstimatorWorkflow.push_back(&OP_SEQ_FILE);
  abundanceEstimatorWorkflow.push_back(&OP_COUNT_FILE);
  abundanceEstimatorWorkflow.push_back(&OP_OUTPREFIX);
  abundanceEstimatorWorkflow.push_back(&OP_SPACED_KMER_PATTERN);
  abundanceEstimatorWorkflow.push_back(&OP_SKIP);
  abundanceEstimatorWorkflow.push_back(&OP_COUNT_MODE);
  abundanceEstimatorWorkflow.push_back(&OP_VERBOSE);

  //profile
  profileWorkflow.push_back(&OP_SEQ_FILE);
  profileWorkflow.push_back(&OP_COUNT_FILE);
  profileWorkflow.push_back(&OP_OUTPREFIX);
  profileWorkflow.push_back(&OP_SPACED_KMER_PATTERN);
  profileWorkflow.push_back(&OP_SKIP);
  profileWorkflow.push_back(&OP_COUNT_MODE);
  profileWorkflow.push_back(&OP_VERBOSE);

  //counts2flat
  counts2flatWorkflow.push_back(&OP_COUNT_FILE);
  counts2flatWorkflow.push_back(&OP_COUNT_MODE);
  counts2flatWorkflow.push_back(&OP_OUTPREFIX);
  counts2flatWorkflow.push_back(&OP_SPACED_KMER_PATTERN);

  }

void Options::setDefaults() {

  spacedKmerPattern="11110111111011011101010111011011111101111";
  skip = 10;

  threshold = 3;
  tolerance = 0.01;
  maxTrimLen = 0;
  maxCorrNum = 0;
  updateLookup = false;

  dropLevel1 = 0.33;
  dropLevel2 = 0.33;

  aligned = false;
  softFilter = false;

  countMode = COUNT_MODE_SUM;

  dryRun = false;
  threads = 1; //TODO
  verbose = Info::INFO;
}

void printToolUsage(const Command &command, const int FLAG) {

  std::stringstream usage;
  usage << "coco " << command.cmd << "\n";
  usage << command.descriptLong << "\n\n";
  usage << "© " << command.author << "\n\n";

  usage << "Usage: coco " << command.cmd << " " << command.usage << "\n\n";

  if (FLAG == EXTENDED) {
    const std::vector<cocoOption*> &options = *command.opt;

    size_t maxParamWidth = 0;
    for (size_t idx = 0; idx < options.size(); idx++) {
      maxParamWidth = std::max(strlen(options[idx]->display),maxParamWidth);
    }

    size_t frontParamWidth = 2;
    maxParamWidth+=6;
      size_t descriptionStart = 0;
      std::string paramString;
    for (size_t idx = 0; idx < options.size(); idx++) {
        paramString.clear();
      paramString += std::string(frontParamWidth, ' ') + options[idx]->display + std::string(maxParamWidth < strlen(options[idx]->display)? 1 : maxParamWidth- strlen(options[idx]->display), ' ');

      descriptionStart = paramString.length();
      char *descr = strdup(options[idx]->description);
      char *ptr = strtok(descr, "\n");
      paramString += std::string(ptr);
      ptr = strtok(NULL, "\n");
      while(ptr != NULL) {
          paramString += "\n" + std::string(descriptionStart,' ') + std::string(ptr);
          ptr = strtok(NULL, "\n");
      }


          usage << paramString << "\n";
        //usage << options[idx]->display << "\t" << options[idx]->description << "\n";
    }
  }
  std::cerr << usage.str() << "\n";
}

void Options::parseOptions(int argc, const char *argv[],
                           const Command &command) {
  std::vector<cocoOption*> &options = *command.opt;
  int opt, longIndex = 0;
  extern char *optarg;

  static const struct option longOpts[] = {
    {"help",      no_argument,       NULL, 'h'},
    {"seqfile",   required_argument, NULL, 0},
    {"counts",    required_argument, NULL, 0},
    {"outprefix",    required_argument, NULL, 0},
    {"spaced-pattern",    required_argument, NULL, 0},
    {"skip",    required_argument, NULL, 0},
    {"threshold", required_argument, NULL, 0},
    {"tolerance", required_argument, NULL, 0},
    {"max-corr-num", required_argument, NULL, 0},
    {"max-trim-len", required_argument, NULL, 0},
    {"dry-run",    required_argument, NULL, 0},
    {"drop-level1",    required_argument, NULL, 0},
    {"drop-level2",    required_argument, NULL, 0},
    {"aligned",    no_argument, NULL, 0},
    {"soft",    no_argument, NULL, 0},
    {"count-mode",   required_argument, NULL, 0},
    {"threads",   required_argument, NULL, 0},
    {"verbose",   required_argument, NULL, 0},
    {NULL,        no_argument,       NULL, 0}
  };

  if (argc < 2) {
    printToolUsage(command, SIMPLE);
    EXIT(EXIT_FAILURE);
  }

  while ((opt = getopt_long(argc, (char **) argv, "h", longOpts, &longIndex)) != -1) {
    switch (opt) {
      case 'h':
        printToolUsage(command, EXTENDED);
        EXIT(EXIT_SUCCESS);
      case 0: {
        std::string optname(longOpts[longIndex].name);
        for (size_t idx = 0; idx < options.size(); idx++) {
          if (optname.compare("--help") == 0) {
            printToolUsage(command, EXTENDED);
            EXIT(EXIT_SUCCESS);
          }

          if (optname.compare(options[idx]->name) == 0) {
            if (options[idx]->isSet) {
              Info(Info::ERROR) <<  "ERROR: Duplicate option " << options[idx]->display << "\n\n";
              printToolUsage(command, SIMPLE);
              EXIT(EXIT_FAILURE);
            }

            if (typeid(std::string) == options[idx]->type) {
              if (optarg[0] == '-') {
                Info(Info::ERROR) << "ERROR: Invalid argument " << optarg << " for option " << options[idx]->display << "\n";
                EXIT(EXIT_FAILURE);
              }

              std::string val(optarg);
              if (val.length() != 0) {
                std::string *currVal = ((std::string *) options[idx]->value);
                currVal->assign(val);
                options[idx]->isSet = true;
              }
            } else if (typeid(int) == options[idx]->type) {
              int val = atoi(optarg);
              *((int *) options[idx]->value) = val;
              options[idx]->isSet = true;

            } else if (typeid(unsigned int) == options[idx]->type) {
              int val = atoi(optarg);
              *((int *) options[idx]->value) = val;
              options[idx]->isSet = true;

            } else if (typeid(bool) == options[idx]->type) {

              *((bool *) options[idx]->value) = true;
              options[idx]->isSet = true;

            } else if (typeid(float) == options[idx]->type) {

              *((float *) options[idx]->value) = std::stof(optarg);
              options[idx]->isSet = true;

            } else if (typeid(double) == options[idx]->type) {

              *((double *) options[idx]->value) = std::stod(optarg);
              options[idx]->isSet = true;

            }
            else {
              Info(Info::ERROR) << "ERROR: Wrong option type in parseOptions used for " << options[idx]->display << "\n" \
                                   "Please send an error report to the developers.\n";
              EXIT(EXIT_FAILURE);
            }
            continue;
          }
        }
        break;
      }

      case '?':
        /* error message already printed by getopt function*/
        std::cerr << "\n";
        printToolUsage(command, SIMPLE);
        EXIT(EXIT_FAILURE);
      default:
        abort();
    }
  }

  if (optind != argc) {
    Info(Info::ERROR) << "ERROR: Superfluous Argument " << argv[optind] << "\n\n";
    printToolUsage(command, SIMPLE);
    EXIT(EXIT_FAILURE);
  }

  for (cocoOption *option: options) {
    //Check if option is required for current tool
    if (command.id & option->required) {
      //If required and not set -> Error
      if (!option->isSet) {
        Info(Info::ERROR) << "ERROR: Option " << option->display << " is required for tool " << command.cmd << "\n";
        EXIT(EXIT_FAILURE);
      }
    }
  }

  Info::setVerboseLevel(verbose);

}
