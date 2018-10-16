/* Option handling for COCO tools/workflows */

#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstddef>
#include <vector>
#include <string>
#include <typeinfo>
#include "ToolInfo.h"

enum coco_tools{ABUNDANCE_ESTIMATOR=1, PCREADS=2, COUNTPROFILE=3};

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
         void *val, unsigned long required):
         name(n), uniqid(uid), display(dp), description(descr), type(ty),value(val),isSet(false),required(required){}
};


/* general class <Options>, provides vector of <cocoOption>s for different
 * COCO workflows */
class Options
{
public:

  void setDefaults();
  static Options& getInstance()
  {
    if (instance == NULL)
    {
      initInstance();
    }
    return *instance;
  }
  static void initInstance()
  {
      new Options;
  }

  void parseOptions(int argc, const char* argv[],
                    const ToolInfo& tool);
                       /*size_t requiredParameterCount,
                       bool printParameters = true,
                       int parseFlags = 0,
                       int outputFlags = 0);*/

  std::string seqFile;
  std::string kcFile;
  unsigned int readAvgLen;
  unsigned int kmerWeight;
  unsigned int threads;

  std::vector<cocoOption> empty;
  std::vector<cocoOption> abundanceEstimatorWorkflow;
  std::vector<cocoOption> pcreadsWorkflow;
  std::vector<cocoOption> profileWorkflow;

  OPTION(OP_SEQ_FILE)
  OPTION(OP_KC_FILE)
  OPTION(OP_AVERAGE_LENGTH)
  OPTION(OP_KMER_WEIGHT)
  OPTION(OP_THREADS)

protected:
    Options();
    static Options* instance;
    virtual ~Options() {};

private:
    Options(Options const&);
    void operator=(Options const&);
};

#endif // OPTIONS_H
