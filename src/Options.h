// Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>

#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstddef>
#include <vector>
#include <string>
#include <typeinfo>
#include "Command.h"

struct cocoOption;
#define OPTION(x) const static int x##_ID = __COUNTER__; \
             cocoOption x;

struct cocoOption {
  const char *name;
  const int uniqid;
  const char *display;
  const char *description;
  const std::type_info &type;
  void *value;
  bool isSet;
  bool isFile;
  //unsigned long required;

  cocoOption(int uid, const char *n, const char *dp,
             const char *descr, const std::type_info &ty,
             void *val, bool isFile) :
    name(n), uniqid(uid), display(dp), description(descr), type(ty), value(val), isSet(false), isFile(isFile) {}
};


class Options {
public:

  void setDefaults();

  static Options &getInstance() {
    if (instance == NULL) {
      initInstance();
    }
    return *instance;
  }

  static void initInstance() {
    new Options;
  }

  static void deleteInstance() {
    delete instance;
  }

  void parseOptions(int argc, const char *argv[], const Command &command);

  void printParameterSettings(const Command &command);

  /*size_t requiredParameterCount,
  bool printParameters = true,
  int parseFlags = 0,
  int outputFlags = 0);*/

  std::string reads;
  std::string forwardReads;
  std::string reverseReads;
  std::string countFile;
  std::string outdir;
  std::string outprefix;

  std::string spacedKmerPattern;
  //int kmerWeight;

  int countMode;
  int skip;

  double threshold;
  int pseudocount;
  int lowerBound;
  int maxTrimLen;
  int maxCorrNum;
  bool updateLookup;

  double dropLevel1;
  double dropLevel2;

  int threads;
  int verbose;

  bool aligned;
  bool softFilter;

  std::vector<cocoOption*> empty;
  std::vector<cocoOption*> filterWorkflow;
  std::vector<cocoOption*> profileWorkflow;
  std::vector<cocoOption*> abundanceEstimatorWorkflow;
  std::vector<cocoOption*> correctionWorkflow;
  std::vector<cocoOption*> consensusWorkflow;
  std::vector<cocoOption*> counts2flatWorkflow;

  OPTION(OP_READS)
  OPTION(OP_FORWARD_READS)
  OPTION(OP_REVERSE_READS)
  OPTION(OP_COUNT_FILE)
  OPTION(OP_OUTDIR)
  OPTION(OP_OUTPREFIX)
  OPTION(OP_SPACED_KMER_PATTERN)
  OPTION(OP_SKIP)
  OPTION(OP_THRESHOLD)
  OPTION(OP_PSEUDOCOUNT)
  OPTION(OP_LOWER_BOUND)
  OPTION(OP_MAX_CORR_NUM)
  OPTION(OP_MAX_TRIM_LEN)
  OPTION(OP_UPDATE_LOOKUPTABLE)
  OPTION(OP_DROP_LEVEL1)
  OPTION(OP_DROP_LEVEL2)
  OPTION(OP_SOFT_FILTER)
  OPTION(OP_ALIGNED)
  OPTION(OP_THREADS)
  OPTION(OP_VERBOSE)
  OPTION(OP_COUNT_MODE)

  static const int COUNT_MODE_SUM = 0;
  static const int COUNT_MODE_MAX = 1;

protected:
  Options();

  static Options *instance;

  virtual ~Options() {};

private:
  Options(Options const &);

  void operator=(Options const &);
};

#endif // OPTIONS_H
