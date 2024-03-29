// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>
#include <iostream>
#include <sstream>
#include <climits>
#include <cstring>
#include <libgen.h>

#include "Options.h"
#include "Info.h"
#include "util.h"
#include "filehandling.h"


Options *Options::instance = NULL;

Options::Options() :
  OP_READS(OP_READS_ID, "--reads", "Unpaired Reads",
               "file with unpaired/single/merged reads (fasta/fastq format)",
               typeid(std::string), (void *) &reads, 1),
  OP_FORWARD_READS(OP_FORWARD_READS_ID, "-1", "Forward Reads",
              "file with forward paired-end reads (fasta/fastq format)",
              typeid(std::string), (void *) &forwardReads, 1),
  OP_REVERSE_READS(OP_REVERSE_READS_ID, "-2", "Reverse Reads",
              "file with reverse paired-end reads (fasta/fastq format)",
              typeid(std::string), (void *) &reverseReads, 1),
  OP_COUNT_FILE(OP_COUNT_FILE_ID, "--counts", "Count File",
               "pre computed solid kmer count file in hdf5 format (dsk output)",
               typeid(std::string), (void *) &countFile, 1),
  OP_OUTDIR(OP_OUTDIR_ID, "--outdir", "Outdir",
               "output directory",
               typeid(std::string), (void *) &outdir, 0),
  OP_OUTPREFIX(OP_OUTPREFIX_ID, "--outprefix", "Outprefix",
               "prefix to use for resultfile(s)",
               typeid(std::string), (void *) &outprefix, 0),
  OP_SPACED_KMER_PATTERN(OP_SPACED_KMER_PATTERN_ID,"--spaced-pattern",
               "Spaced pattern", "user-specified spaced k-mer pattern\n(span must be <=64, 12 <= weight <=32 and symmetric)",
               typeid(std::string), (void*) &spacedKmerPattern, 0),
  OP_SKIP(OP_SKIP_ID,"--skip", "skip", "skip sequences with less than this many k-mers",
               typeid(int), (void*) &skip, 0),
  OP_THRESHOLD(OP_THRESHOLD_ID, "--threshold", "Threshold",
               "percent drop threshold relative to neighborhood counts",
               typeid(double), (void *) &threshold, 0),
  OP_PSEUDOCOUNT(OP_PSEUDOCOUNT_ID, "--pseudocount", "Pseudocount",
               "untrusted count added to the pseudocount parameter",
               typeid(int), (void *) &pseudocount, 0),
  OP_LOWER_BOUND(OP_LOWER_BOUND_ID, "--lowerbound", "Lower bound",
                 "lower bound for neighborhood counts to be considered",
                 typeid(int), (void *) &lowerBound, 0),
  OP_MAX_CORR_NUM(OP_MAX_CORR_NUM_ID, "--max-corr-num", "Max number of correction per read",
               "maximal number of corrections performed per read, changes are discarded otherwise",
               typeid(int), (void *) &maxCorrNum, 0),
  OP_MAX_TRIM_LEN(OP_MAX_TRIM_LEN_ID, "--max-trim-len", "Max number of trimmed nucleotides",
               "maximal number of nucleotides trimmed from the beginning or end of a read if error could not be corrected",
               typeid(int), (void *) &maxTrimLen, 0),
  OP_UPDATE_LOOKUPTABLE(OP_UPDATE_LOOKUPTABLE_ID, "--update-lookup", "Update lookup table",
               "update counts in lookuptable after a sequence is corrected\n"\
               "(slow down and creates read order dependency but might help in low coverage regions)",
                typeid(bool), (void *) &updateLookup, 0),
  OP_DROP_LEVEL1(OP_DROP_LEVEL1_ID, "--drop-level1", "Drop-level1",
               "local drop criterion (range 0-0.33)",
               typeid(double), (void *) &dropLevel1, 0),
  OP_DROP_LEVEL2(OP_DROP_LEVEL2_ID, "--drop-level2", "Drop-level2",
               "global drop criterion (range 0-0.33)",
               typeid(double), (void *) &dropLevel2, 0),
  OP_SOFT_FILTER(OP_SOFT_FILTER_ID, "--soft", "Soft Filtering",
               "less strict filtering mode due to more strict masking strategy ",
               typeid(bool), (void *) &softFilter, 0),
  OP_ALIGNED(OP_ALIGNED_ID, "--aligned", "Aligned",
               "optimize abundance estimation for reads that span the same region (amplicon sequence data)",
               typeid(bool), (void *) &aligned, 0),
  OP_THREADS(OP_THREADS_ID, "--threads", "Number of threads", "number of threads, not supported yet", typeid(int), (void *) &threads, 0),
  OP_VERBOSE(OP_VERBOSE_ID, "--verbose", "Verbosity level", "verbosity level, 0: quiet 1: Errors, 2: +Warnings, 3: +Info, 4: +Debug",
             typeid(int), (void *) &verbose, 0),
  // expert options
  OP_COUNT_MODE(OP_COUNT_MODE_ID, "--count-mode", "Count mode",
                "way to store counts for concurrent spaced kmers (expert option)\n 0: sum \n 1: maximize",
                typeid(int), (void *) &countMode, 0)

  {
  if (instance) {
    std::cerr << "Parameter instance already exists!\n";
    abort();
  }
  instance = this;

  setDefaults();

  //TODO: threads

  //correction
  correctionWorkflow.push_back(&OP_FORWARD_READS);
  correctionWorkflow.push_back(&OP_REVERSE_READS);
  correctionWorkflow.push_back(&OP_READS);
  correctionWorkflow.push_back(&OP_COUNT_FILE);
  correctionWorkflow.push_back(&OP_COUNT_MODE);
  correctionWorkflow.push_back(&OP_OUTDIR);
  correctionWorkflow.push_back(&OP_OUTPREFIX);
  correctionWorkflow.push_back(&OP_SPACED_KMER_PATTERN);
  correctionWorkflow.push_back(&OP_SKIP);
  correctionWorkflow.push_back(&OP_THRESHOLD);
  correctionWorkflow.push_back(&OP_PSEUDOCOUNT);
  correctionWorkflow.push_back(&OP_LOWER_BOUND);
  correctionWorkflow.push_back(&OP_MAX_CORR_NUM);
  correctionWorkflow.push_back(&OP_MAX_TRIM_LEN);
  correctionWorkflow.push_back(&OP_UPDATE_LOOKUPTABLE);
  correctionWorkflow.push_back(&OP_VERBOSE);

  // filter
  filterWorkflow.push_back(&OP_FORWARD_READS);
  filterWorkflow.push_back(&OP_REVERSE_READS);
  filterWorkflow.push_back(&OP_READS);
  filterWorkflow.push_back(&OP_COUNT_FILE);
  filterWorkflow.push_back(&OP_COUNT_MODE);
  filterWorkflow.push_back(&OP_OUTDIR);
  filterWorkflow.push_back(&OP_OUTPREFIX);
  filterWorkflow.push_back(&OP_SPACED_KMER_PATTERN);
  filterWorkflow.push_back(&OP_SKIP);
  filterWorkflow.push_back(&OP_THRESHOLD);
  filterWorkflow.push_back(&OP_VERBOSE);

  //abundanceEstimator
  abundanceEstimatorWorkflow.push_back(&OP_FORWARD_READS);
  abundanceEstimatorWorkflow.push_back(&OP_REVERSE_READS);
  abundanceEstimatorWorkflow.push_back(&OP_READS);
  abundanceEstimatorWorkflow.push_back(&OP_COUNT_FILE);
  abundanceEstimatorWorkflow.push_back(&OP_COUNT_MODE);
  abundanceEstimatorWorkflow.push_back(&OP_OUTDIR);
  abundanceEstimatorWorkflow.push_back(&OP_OUTPREFIX);
  abundanceEstimatorWorkflow.push_back(&OP_SPACED_KMER_PATTERN);
  abundanceEstimatorWorkflow.push_back(&OP_SKIP);
  abundanceEstimatorWorkflow.push_back(&OP_VERBOSE);

  //profile
  profileWorkflow.push_back(&OP_FORWARD_READS);
  profileWorkflow.push_back(&OP_REVERSE_READS);
  profileWorkflow.push_back(&OP_READS);
  profileWorkflow.push_back(&OP_COUNT_FILE);
  profileWorkflow.push_back(&OP_COUNT_MODE);
  profileWorkflow.push_back(&OP_OUTDIR);
  profileWorkflow.push_back(&OP_OUTPREFIX);
  profileWorkflow.push_back(&OP_SPACED_KMER_PATTERN);
  profileWorkflow.push_back(&OP_SKIP);
  profileWorkflow.push_back(&OP_VERBOSE);

  //counts2flat
  counts2flatWorkflow.push_back(&OP_COUNT_FILE);
  counts2flatWorkflow.push_back(&OP_COUNT_MODE);
  counts2flatWorkflow.push_back(&OP_OUTDIR);
  counts2flatWorkflow.push_back(&OP_OUTPREFIX);
  counts2flatWorkflow.push_back(&OP_SPACED_KMER_PATTERN);

  }

void Options::setDefaults() {

  reads = "";
  forwardReads = "";
  reverseReads = "";

  outprefix="";
  outdir="coco_out/";

  spacedKmerPattern="11110111111011011101010111011011111101111";
  skip = 10;

  pseudocount = 1;
  threshold = 0.01;
  lowerBound = 5;
  maxTrimLen = 0;
  maxCorrNum = 10;
  updateLookup = false;

  dropLevel1 = 0.33;
  dropLevel2 = 0.33;

  aligned = false;
  softFilter = false;

  countMode = COUNT_MODE_SUM;

  threads = 1; //TODO
  verbose = Info::INFO;
}

void printToolUsage(const Command &command, const int FLAG) {

  std::stringstream usage;
  /*usage << "coco " << command.cmd << "\n";
  usage << command.descriptLong << "\n\n";*/
  usage << "© " << command.author << "\n\n";

  usage << "Usage: coco " << command.cmd << " " << command.usage << "\n\n";

  if (FLAG == EXTENDED) {
    const std::vector<cocoOption*> &options = *command.opt;

    size_t maxParamWidth = 0;
    for (size_t idx = 0; idx < options.size(); idx++) {
      maxParamWidth = std::max(strlen(options[idx]->name),maxParamWidth);
    }

    size_t frontParamWidth = 2;
    maxParamWidth+=6;
    std::string paramString;
    for (size_t idx = 0; idx < options.size(); idx++) {
        paramString.clear();
      paramString += std::string(frontParamWidth, ' ') + options[idx]->name +
                     std::string(maxParamWidth < strlen(options[idx]->name)? 1 : maxParamWidth- strlen(options[idx]->name), ' ');

      size_t descriptionStart = paramString.length();
      char *descr = strdup(options[idx]->description);
      char *ptr = strtok(descr, "\n");
      paramString += std::string(ptr);

      std::string valueString;
      if (typeid(std::string) == options[idx]->type) {
        valueString = *((std::string *) options[idx]->value);
      } else if (typeid(int) == options[idx]->type) {
        valueString = std::to_string(*(int *) options[idx]->value);
      } else if (typeid(unsigned int) == options[idx]->type) {
        valueString = std::to_string(*(int *) options[idx]->value);
      } else if (typeid(bool) == options[idx]->type) {
        valueString = std::to_string(*(bool *) options[idx]->value);
      } else if (typeid(float) == options[idx]->type) {
        //valueString = std::to_string(*(float *) options[idx]->value);
        char buffer[32];
        int n = sprintf(buffer, "%.3f", *(float *) options[idx]->value);
        valueString = std::string(buffer, n);
      } else if (typeid(double) == options[idx]->type) {
        //valueString = std::to_string(*(double *) options[idx]->value);
        char buffer[32];
        int n = sprintf(buffer, "%.3lf", *(double *) options[idx]->value);
        valueString = std::string(buffer, n);
      }
      if (valueString.length() > 0) {
        paramString += std::string(" [") + valueString + std::string("]");
      }
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

void Options::printParameterSettings(const Command &command) {

  std::stringstream settings;
  const std::vector<cocoOption*> &options = *command.opt;

  size_t maxParamWidth = 0;
  for (size_t idx = 0; idx < options.size(); idx++) {
    maxParamWidth = std::max(strlen(options[idx]->display),maxParamWidth);
  }

  size_t frontParamWidth = 2;
  maxParamWidth+=6;
  std::string paramString;
  for (size_t idx = 0; idx < options.size(); idx++) {
      paramString.clear();
      paramString += std::string(frontParamWidth, ' ') + options[idx]->display +
                     std::string(maxParamWidth < strlen(options[idx]->display)? 1 : maxParamWidth- strlen(options[idx]->display), ' ');

      std::string valueString;
      if (typeid(std::string) == options[idx]->type) {
        valueString = *((std::string *) options[idx]->value);
      } else if (typeid(int) == options[idx]->type) {
        valueString = std::to_string(*(int *) options[idx]->value);
      } else if (typeid(unsigned int) == options[idx]->type) {
        valueString = std::to_string(*(int *) options[idx]->value);
      } else if (typeid(bool) == options[idx]->type) {
        valueString = std::to_string(*(bool *) options[idx]->value);
      } else if (typeid(float) == options[idx]->type) {
        valueString = std::to_string(*(float *) options[idx]->value);
      } else if (typeid(double) == options[idx]->type) {
        valueString = std::to_string(*(double *) options[idx]->value);
      }
      if (valueString.length() > 0) {
        //paramString += valueString;
        paramString += std::string(" [") + valueString + std::string("]");
      }

    settings << paramString << "\n";
  }
  Info(Info::INFO) << settings.str() << "\n";
}

void Options::parseOptions(int argc, const char *argv[], const Command &command) {
  std::vector<cocoOption*> &options = *command.opt;

  if (argc < 1) {
    printToolUsage(command, EXTENDED);
    EXIT(EXIT_FAILURE);
  }

 // if -h or --help is set print usage without parsing other arguments
  for (int argIdx = 0; argIdx < argc; argIdx++) {

    std::string optname(argv[argIdx]);
    if ((optname.compare("-h") == 0 || optname.compare("--help") == 0)) {
      printToolUsage(command, EXTENDED);
      EXIT(EXIT_SUCCESS);
    }
  }

  // parse other arguments
  for (int argIdx = 0; argIdx < argc; argIdx++) {

    std::string optname(argv[argIdx]);
    bool recognizedParameter = false;

    for (size_t idx = 0; idx < options.size(); idx++) {

      if (optname.compare(options[idx]->name) == 0) {
        recognizedParameter = true;
        if (options[idx]->isSet) {
          Info(Info::ERROR) <<  "ERROR: Duplicate option " << options[idx]->name << "\n\n";
          printToolUsage(command, SIMPLE);
          EXIT(EXIT_FAILURE);
        }

        if (typeid(bool) != options[idx]->type && (argIdx == argc-1 || argv[argIdx+1][0] == '-')) {
          Info(Info::ERROR) << "ERROR: No argument for option " << options[idx]->name << " is given \n";
          EXIT(EXIT_FAILURE);
        }
        const char *optarg = argv[argIdx+1];
        if (typeid(std::string) == options[idx]->type) {

          std::string val(optarg);
          std::string *currVal = ((std::string *) options[idx]->value);
          currVal->assign(val);
          options[idx]->isSet = true;
          argIdx++;

        } else if (typeid(int) == options[idx]->type) {
          int val = atoi(optarg);
          *((int *) options[idx]->value) = val;
          options[idx]->isSet = true;
          argIdx++;

        } else if (typeid(unsigned int) == options[idx]->type) {
          int val = atoi(optarg);
          *((int *) options[idx]->value) = val;
          options[idx]->isSet = true;
          argIdx++;

        } else if (typeid(bool) == options[idx]->type) {

          if(argv[argIdx+1][0] == '-'){ //no argument for bool option -> true
            *((bool *) options[idx]->value) = true;
          }else {
            if (strcmp(optarg, "true") == 0 || strcmp(optarg, "TRUE") == 0 || strcmp(optarg, "1") == 0) {
              *((bool *) options[idx]->value) = true;
            } else if (strcmp(optarg, "false") == 0 || strcmp(optarg, "FALSE") == 0 || strcmp(optarg, "0") == 0) {
              *((bool *) options[idx]->value) = false;
            } else {
              Info(Info::ERROR) << "ERROR: Invalid boolean argument for option " << options[idx]->name << "\n";
              EXIT(EXIT_FAILURE);
            }
            argIdx++;
          }
          options[idx]->isSet = true;
        } else if (typeid(float) == options[idx]->type) {

          *((float *) options[idx]->value) = std::stof(optarg);
          options[idx]->isSet = true;
          argIdx++;

        } else if (typeid(double) == options[idx]->type) {

          *((double *) options[idx]->value) = std::stod(optarg);
          options[idx]->isSet = true;
          argIdx++;

        }
        else {
          Info(Info::ERROR) << "ERROR: Wrong option type in parseOptions used for " << options[idx]->name << "\n" \
                               "Please send an error report to the developers.\n";
          EXIT(EXIT_FAILURE);
        }
        break;
      }
    }
    if (!recognizedParameter){
      Info(Info::ERROR) << "ERROR: Unrecognized parameter " << optname << "\n";
      printToolUsage(command, EXTENDED);
      EXIT(EXIT_FAILURE);
    }
  }

  // required for all tools except counts2flat
  if(command.id!=COUNTS2FLAT && (!((this->OP_READS.isSet)|(this->OP_FORWARD_READS.isSet && this->OP_REVERSE_READS.isSet)))){
    Info(Info::ERROR) << "ERROR: Either " << this->OP_READS.name << " or " << this->OP_FORWARD_READS.name << " and "
                      << this->OP_REVERSE_READS.name << " must be set\n";
    EXIT(EXIT_FAILURE);
  } else if (command.id==COUNTS2FLAT) {
    if(!this->OP_COUNT_FILE.isSet) {
      Info(Info::ERROR) << "ERROR: " << this->OP_COUNT_FILE.name << " must be set\n";
      EXIT(EXIT_FAILURE);
    }
  }

  for (cocoOption *option: options) {
    if (option->isSet && option->isFile) {
      if (!fileExists(((std::string *) option->value)->c_str())) {
        Info(Info::ERROR) << "ERROR: In Option " << option->name << " file '" \
                          << *(std::string *) option->value << "' does not exist\n";
        EXIT(EXIT_FAILURE);
      }
    }
  }

  if(outdir.length() > 0) {
    if (outdir[outdir.length() - 1] != '/')
      outdir += "/";

    if (directoryExists(outdir.c_str()) == false) {
      if (_mkdir(outdir) == false) {
        Info(Info::ERROR) << "ERROR: Failed to create output directory " << outdir << "\n";
        EXIT(EXIT_FAILURE);
      }
    }
  }

  if (this->OP_OUTPREFIX.isSet) {
    char tmp[1024];
    strcpy(tmp, outprefix.c_str());
    if (strcmp(dirname(tmp), ".") != 0) {
      Info(Info::ERROR) << "ERROR: Value for option " << this->OP_OUTPREFIX.name << " is a path, please use "
                        << this->OP_OUTDIR.name << " to set a directory\n";
      EXIT(EXIT_FAILURE);
    }
  }



  /*for (cocoOption *option: options) {
    //Check if option is required for current tool
    if (command.id & option->required) {
      //If required and not set -> Error
      if (!option->isSet) {
        Info(Info::ERROR) << "ERROR: Option " << option->name << " is required for " << command.cmd << " command, "\
        "but not set\n";
        EXIT(EXIT_FAILURE);
      }
    }
  }*/

  Info::setVerboseLevel(verbose);

}
