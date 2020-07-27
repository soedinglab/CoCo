// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include "options.h"
#include "util.h"

Options *Options::instance = NULL;

Options::Options() :
  OP_SEQ_FILE(OP_SEQ_FILE_ID, "seqFile", "--seqFile",
              "sequence file (reads or contigs in fasta format)",
              typeid(std::string), (void *) &seqFile, PROFILE | FILTER | ABUNDANCE_ESTIMATOR | CONSENSUS),
  OP_COUNT_FILE(OP_COUNT_FILE_ID, "counts", "--counts ",
                "pre computed kmer count file in hdf5 format (dsk output format), Note: only supports 41-mers yet",
                typeid(std::string), (void *) &countFile, 0),
  OP_OUTPREFIX(OP_OUTPREFIX_ID, "outprefix", "--outprefix ",
               "prefix to use for resultfile(s)",
               typeid(std::string), (void *) &outprefix, 0),
  /*
  OP_KMER_WEIGHT(OP_KMER_WEIGHT_ID,"kmerWeight",
  "--kmerWeight", "number of informative positions in a k-mer pattern, "
  "default: 27", typeid(int),  (void *) &kmerWeight,0)
  ,*/
  OP_THREADS(OP_THREADS_ID, "threads", "--threads", "number of threads"
                                                    "default: 1", typeid(int), (void *) &threads, 0) {
  if (instance) {
    std::cerr << "Parameter instance already exists!\n";
    abort();
  }
  instance = this;

  setDefaults();

  //profile
  profileWorkflow.push_back(OP_SEQ_FILE);
  profileWorkflow.push_back(OP_COUNT_FILE);
  profileWorkflow.push_back(OP_THREADS);
  profileWorkflow.push_back(OP_OUTPREFIX);

  // filter
  filterWorkflow.push_back(OP_SEQ_FILE);
  filterWorkflow.push_back(OP_COUNT_FILE);
  filterWorkflow.push_back(OP_THREADS);
  filterWorkflow.push_back(OP_OUTPREFIX);

  //abundanceEstimator
  abundanceEstimatorWorkflow.push_back(OP_SEQ_FILE);
  abundanceEstimatorWorkflow.push_back(OP_COUNT_FILE);
  abundanceEstimatorWorkflow.push_back(OP_THREADS);
  abundanceEstimatorWorkflow.push_back(OP_OUTPREFIX);

}

void Options::setDefaults() {
  threads = 1; //TODO
}

void printToolUsage(const Command &command, const int FLAG) {
  std::stringstream usage;
  usage << "coco " << command.cmd << "\n";
  usage << command.descriptLong << "\n\n";
  usage << "© " << command.author << "\n\n";

  usage << "Usage: coco " << command.cmd << " " << command.usage << "\n\n";

  if (FLAG == EXTENDED) {
    const std::vector<cocoOption> &options = *command.opt;
    for (size_t idx = 0; idx < options.size(); idx++) {
      //TODO: layout, maxwidth + 2 space in the beginning
      usage << options[idx].display << "\t" << options[idx].description << "\n";
    }
  }
  std::cerr << usage.str() << "\n";
}

void Options::parseOptions(int argc, const char *argv[],
                           const Command &command) {
  std::vector<cocoOption> &options = *command.opt;
  int opt, longIndex = 0;
  extern char *optarg;

  static const struct option longOpts[] = {
    {"help",      no_argument,       NULL, 'h'},
    {"seqFile",   required_argument, NULL, 0},
    {"countFile", required_argument, NULL, 0},
    {"threads",   required_argument, NULL, 0},
    {NULL,        no_argument,       NULL, 0}
  };

  while ((opt = getopt_long(argc, (char **) argv, "h", longOpts, &longIndex)) != -1) {
    switch (opt) {
      //fprintf(stderr, "longidx=%d", longIndex);
      //fprintf(stderr, "processing %s\n", longOpts[longIndex].name);
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

          if (optname.compare(options[idx].name) == 0) {
            if (options[idx].isSet) {
              fprintf(stderr, "Duplicate option %s\n\n", options[idx].display);
              printToolUsage(command, SIMPLE);
              EXIT(EXIT_FAILURE);
            }

            if (typeid(std::string) == options[idx].type) {
              if (optarg[0] == '-') {
                fprintf(stderr, "Invalid argument %s for option %s\n",
                        optarg, options[idx].display);
                EXIT(EXIT_FAILURE);
              }
              std::string val(optarg);
              if (val.length() != 0) {
                std::string *currVal = ((std::string *) options[idx].value);
                currVal->assign(val);
                options[idx].isSet = true;
              }
            } else if (typeid(int) == options[idx].type) {
              int val = atoi(optarg);
              if (val == 0) {
                fprintf(stderr, "Invalid argument %s for option %s\n",
                        optarg, options[idx].display);
                EXIT(EXIT_FAILURE);
              }
              *((int *) options[idx].value) = val;
              options[idx].isSet = true;

            } else {
              fprintf(stderr, "Wrong option type in parseOptions. "\
                         "please send an error report to the developers.");
              EXIT(EXIT_FAILURE);
            }
            continue;
          }
        }
        //assert(false); //TODO
        break;
      }

      case '?':
        /* error message already printed by getopt function*/
        printToolUsage(command, SIMPLE);
        EXIT(EXIT_FAILURE);
      default:
        abort();
    }
  }

  if (optind != argc) {
    fprintf(stderr, "ERROR: Superfluous Argument %s\n\n", argv[optind]);
    printToolUsage(command, SIMPLE);
    EXIT(EXIT_FAILURE);
  }

  //TODO: valid range?
  for (cocoOption option: options) {
    //Check if option is required for current tool
    if (command.id & option.required) {
      //If required and not set -> Error
      if (!option.isSet) {
        fprintf(stderr, "ERROR: Option %s is required for %s\n", option.display, command.cmd);
        EXIT(EXIT_FAILURE);
      }
    }
  }
  //TODO: flag für print Parameter
}
