#include <cstdio>
#include "ToolInfo.h"
#include "options.h"

int pcoverage(int argc, const char **argv, const ToolInfo* tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

  printf("argument: sampleListFile: %s\n", opt.sampleListFile);

  return 0;
}
