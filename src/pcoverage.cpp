/*
 * Copyright (C) 2018 Annika Seidel <annika.seidel@mpibpc.mpg.de>
 *
 */

#include <cstdio>
#include "ToolInfo.h"
#include "options.h"
#include "types.h"

int pcoverage(int argc, const char **argv, const ToolInfo* tool)
{
  Options &opt = Options::getInstance();
  opt.parseOptions(argc, argv, *tool);

  printf("argument: sampleList: %s\n", opt.sampleListFile.c_str());
  printf("argument: kmerWeight: %u\n", opt.kmerWeight);

  initialize();
  return 0;
}
