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
#define EXIT(exitCode) do { int __status = (exitCode); std::cerr.flush(); std::cout.flush(); exit(__status); } while(0)
#endif

const char* binary_name = "coco";
const char* tool_name = "COCO";
const char* tool_introduction = "COnsensus COmputation";
const char* main_author = "Annika Seidel (annika.seidel@mpibpc.mpg.de)";

extern int pcreads(int argc, const char **argv, const struct ToolInfo* tool);
extern int pcoverage(int argc, const char **argv, const struct ToolInfo* tool);

//Parameters& par = Parameters::getInstance();
struct ToolInfo tool1 = {"pcoverage", pcoverage, "calculates for every read an estimated value for "\
                                        "the population coverage",
                "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
                "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <i:kmer-countFile1.hdf5>...<i:kmer-countFileN.hdf5>",
               };
struct ToolInfo tool2 = {"pcreads", pcreads, "calculates for every read the consensus read",
                     "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
                     "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <i:kmer-countFile.hdf5>",
                    };

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

int main(int argc, const char * argv[])
{
  //TODO: check 64 system, avx2, sse3

  if (argc < 2)
  {
    //TODO: print usage
    return(EXIT_FAILURE);
  }

  if (argv[1][0] == '-' && argv[1][1] == 'h')
  {
    //TODO: print help
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
    //TODO: invalid command, print usage + help
    return(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
