// Written by Annika Seidel

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <sys/stat.h>

#include <gatb/gatb_core.hpp>
#include "kseq.h"

#include "src/ToolInfo.h"
#ifndef EXIT
#define EXIT(exitCode) do { int __status = (exitCode);\
                            std::cerr.flush(); std::cout.flush();\
                            exit(__status); } while(0)
#endif


#ifdef GIT_SHA1
#define str2(s) #s
#define str(s) str2(s)
    const char *version = str(GIT_SHA1);
#undef str
#undef str2
#else
    const char *version = "UNKNOWN";
#endif


const char* binary_name = "coco";
const char* tool_name = "COCO";
const char* tool_introduction = "COCO is an open-source software suite for "\
                                "different COnsensus COmputation applications "\
                                "on short reads";
const char* main_author = "Annika Seidel (annika.seidel@mpibpc.mpg.de)";

extern int pcreads(int argc, const char **argv, const struct ToolInfo* tool);
extern int pcoverage(int argc, const char **argv, const struct ToolInfo* tool);

//Parameters& par = Parameters::getInstance();

std::vector<struct ToolInfo> tools =
{
  {"pcoverage", pcoverage, "calculates for every read an estimated value for "\
                           "the population coverage",
   "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
   "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <i:kmer-countFile1.hdf5>...<i:kmer-countFileN.hdf5>",
  },
  {"pcreads", pcreads, "calculates for every read the consensus read",
   "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
   "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <i:kmer-countFile.hdf5>"
  },
};

struct ToolInfo *getToolInfo(const char *name)
{
  for(size_t idx=0; idx<tools.size(); idx++)
  {
    struct ToolInfo *tool = &tools[idx];
    if(!strcmp(name, tool->cmd))
      return tool;
  }
  return NULL;
}

void printUsage()
{
  std::stringstream usage;
  usage << tool_introduction << "\n\n";
  usage << tool_name << " Version: " << version << "\n";
  usage << "© " << main_author << "\n\n";

  for (size_t j = 0; j < tools.size(); j++)
  {
    struct ToolInfo &t = tools[j];
    usage << "  " + std::string(t.cmd) << "\t"\
          << t.descriptShort << "\n";
  }

  std::cerr << usage.str() << "\n";
}

int main(int argc, const char * argv[])
{
  //TODO: check 64 system, avx2, sse3

  if (argc < 2)
  {;
    printUsage();
    return(EXIT_FAILURE);
  }

  if (argv[1][0] == '-' && argv[1][1] == 'h')
  {
    printUsage(); //TODO: usage extend for help
    return(EXIT_SUCCESS);
  }

  //TODO: configure
  struct ToolInfo *tool;
  if((tool = getToolInfo(argv[1]))!=NULL)
  {
    EXIT(tool->callerFunction(argc-2, argv+2, tool));
  }
  else
  {
    std::cerr << "\nInvalid Command: " << argv[1] << "\n";
    printUsage();
    return(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
