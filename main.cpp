/*
 * Copyright (C) 2018 Annika Seidel <annika.seidel@mpibpc.mpg.de>
 *
 */

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

#include "src/util.h"
#include "src/options.h"
#include "src/ToolInfo.h"


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
extern int profile(int argc, const char **argv, const struct ToolInfo* tool);

Options& opt = Options::getInstance();
std::vector<struct ToolInfo> tools =
{
  {"pcoverage", pcoverage, &opt.pcoverageWorkflow,"estimates for every sequence "\
   "(reads or contigs) a estimated abundance value in one or multiple samples.\n\n",
   "Calculates for every sequence in a given file of concatenated samples "
   "S={s_1,...,s_n} an estimated value for the abundance in one or multiple samples "\
   "t_1,...,t_m. In general {t1,...,tm,} is a subset of S. Provide for every "
   "sample in S the reads or contigs in a concatenated fasta/fastq format and "
   "for every sample in T a kmer-count File in hdf5 format.\n\n"
   "Further use case: abundance values of reads or contigs through many samples "
   "enable binning steps. Call tool several times with same sequence file but "
   "different h5 sample files or use the tool with the multi h5 option (TODO!!!) "
   "to get abundance values across many samples. Join them with TODO to one big"
   "matrix abundance file and provide this for the binning step",
   "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
   "--sampleList <fileWithSampleNamesofS> --kmerCountList<fileWithSampleNamesofT>",
   PCOVERAGE //tool enum in option.h
  },
  {"pcreads", pcreads, &opt.pcreadsWorkflow, "calculates for every read the "\
   "consensus read", "TODO: long discreption",
   "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
   "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <i:kmer-countFile.hdf5>",
   PCREADS //tool enum in option.h
  },
  {"countprofile", profile, &opt.profileWorkflow, "write for every read the "\
   "the maximized spaced k-mer count profile in a seperated file",
   "TODO: long discreption",
   "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
   "<i:fastaFile1[.gz]> ... <i:fastaFileN[.gz]> <i:kmer-countFile.hdf5>",
   PCREADS //tool enum in option.h
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

void printUsage(const int mode=SIMPLE)
{
  std::stringstream usage;
  usage << tool_introduction << "\n\n";
  usage << tool_name << " Version: " << version << "\n";
  usage << "© " << main_author << "\n\n";

  usage << "available tools:\n";
  for (size_t j = 0; j < tools.size(); j++)
  {
    struct ToolInfo &t = tools[j];
    usage << "  " + std::string(t.cmd) << "\t"\
          << t.descriptShort << "\n";
  }

  //TODO: add EXTENDED mode for special tools, show when -h/--help is used

  std::cerr << usage.str() << "\n";
}

int main(int argc, const char * argv[])
{
  //TODO: check 64 system, avx2, sse3

  if (argc < 2)
  {
    printUsage(SIMPLE);
    return(EXIT_FAILURE);
  }

  if ((argv[1][0] == '-' && argv[1][1] == 'h') ||
      strcmp(argv[1], "--help")==0)
  {
    printUsage(EXTENDED); //TODO: usage extend for help
    return(EXIT_SUCCESS);
  }

  //TODO: configure
  struct ToolInfo *tool;
  if((tool = getToolInfo(argv[1]))!=NULL)
  {
    fprintf(stdout, "execute tool: %s\n", tool->cmd);
    EXIT(tool->callerFunction(argc-1, argv+1, tool));
  }
  else
  {
    fprintf(stderr ,"Invalid Command: %s\n", argv[1]);
    printUsage(SIMPLE);
    return(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
