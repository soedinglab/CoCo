#include <cstdio>
#include "ToolInfo.h"
#include "options.h"

int pcoverage(int argc, const char **argv, const ToolInfo* tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

  printf("argument: sampleList: %s\n", opt.sampleListFile);
  printf("argument: kmerWeight: %u\n", opt.kmerWeight);

  return 0;
}
