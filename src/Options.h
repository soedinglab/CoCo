// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

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
  unsigned long required;

  cocoOption(int uid, const char *n, const char *dp,
             const char *descr, const std::type_info &ty,
             void *val, unsigned long required) :
    name(n), uniqid(uid), display(dp), description(descr), type(ty), value(val), isSet(false), required(required) {}
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

  void parseOptions(int argc, const char *argv[],
                    const Command &command);

  /*size_t requiredParameterCount,
  bool printParameters = true,
  int parseFlags = 0,
  int outputFlags = 0);*/

  std::string seqFile;
  std::string countFile;
  std::string outprefix;
  float dropLevel1;
  float dropLevel2;
  //unsigned int kmerWeight;
  unsigned int threads;
  int stepsize;
  int span;
  int weight;
  unsigned long long int pmstart;
  unsigned long long int pmstop;
  int rand;
  int countMode;
  int verbose;
  bool aligned;
  bool softFilter;
  bool dryRun;

  std::vector<cocoOption*> empty;
  std::vector<cocoOption*> filterWorkflow;
  std::vector<cocoOption*> profileWorkflow;
  std::vector<cocoOption*> abundanceEstimatorWorkflow;
  std::vector<cocoOption*> correctionWorkflow;
  std::vector<cocoOption*> consensusWorkflow;

  OPTION(OP_SEQ_FILE)
  OPTION(OP_COUNT_FILE)
  OPTION(OP_OUTPREFIX)
  OPTION(OP_DRY_RUN)
  OPTION(OP_DROP_LEVEL1)
  OPTION(OP_DROP_LEVEL2)
  OPTION(OP_ALIGNED)
  OPTION(OP_SOFT_FILTER)
  OPTION(OP_COUNT_MODE)
  OPTION(OP_THREADS)
  OPTION(OP_VERBOSE)
  OPTION(OP_STEPSIZE)
  OPTION(OP_SPAN)
  OPTION(OP_WEIGHT)
  OPTION(OP_PMSTART)
  OPTION(OP_PMSTOP)
  OPTION(OP_RAND)

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
