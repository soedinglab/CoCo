/* Option handling for COCO tools/workflows */

#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstddef>
#include <vector>
#include <typeinfo>
#include "ToolInfo.h"

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

  cocoOption(int uid, const char *n, const char *dp,
         const char *descr, const std::type_info &ty,
         void *val):
         uniqid(uid), name(n), display(dp), description(descr), type(ty),value(val),isSet(false){}
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

  const char* sampleListFile;
  const char* kmerCountListFile;
  const char *readAvgLenFile;
  unsigned int kmerWeight;

  std::vector<cocoOption> empty;
  std::vector<cocoOption> pcoverageWorkflow;
  std::vector<cocoOption> pcreadsWorkflow;

  OPTION(OP_SAMPLE_LIST)
  OPTION(OP_KC_LIST)
  OPTION(OP_READ_AVERAGELEN_LIST)
  OPTION(OP_KMER_WEIGHT)

protected:
    Options();
    static Options* instance;
    virtual ~Options() {};

private:
    Options(Options const&);
    void operator=(Options const&);
};

#endif // OPTIONS_H
