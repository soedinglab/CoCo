//Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>

#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>
#include "src/util.h"
#include "src/Options.h"
#include "src/Command.h"
#include "src/Info.h"


#ifdef GIT_SHA1
#define str2(s) #s
#define str(s) str2(s)
  const char *version = str(GIT_SHA1);
#undef str
#undef str2
#else
  const char *version = "UNKNOWN";
#endif

const char *tool_name = "CoCo";
const char *tool_introduction = "CoCo is an open-source software suite for "\
                                "different COnsensus COmputation applications "\
                                "on short reads and contigs";

const char *main_author = "Annika Jochheim (annika.jochheim@mpibpc.mpg.de)";

extern int correction(int argc, const char **argv, const struct Command *tool);

extern int consensus(int argc, const char **argv, const struct Command *tool);

extern int filter(int argc, const char **argv, const struct Command *tool);

extern int abundanceEstimator(int argc, const char **argv, const struct Command *tool);

extern int profile(int argc, const char **argv, const struct Command *tool);

extern int counts2flat(int argc, const char **argv, const struct Command *tool);

Options &opt = Options::getInstance();
std::vector<struct Command> commands =
  {
     /*** basic tools ***/

    {"correction", correction, &opt.correctionWorkflow, "correct sequencing errors",
      "identify reads with sequencing errors and correct them",
      "Annika Jochheim <annika.jochheim@mpibpc.mpg.de>",
      " -1 <fasta|q> -2 <fasta|q> | --reads <fasta|q> [--counts <count.h5>] [options]",
      CORRECTOR
    },
  /*{"consensus", consensus, &opt.consensusWorkflow, "generate consensus reads ",
     "flip SNPs to the major allele to generate the consensus nucleotide sequence",
     "Annika Jochheim <annika.jochheim@mpibpc.mpg.de>",
     "--reads <fastaFile> [--counts <count.h5>] [--outprefix <string>] [options]",
     CONSENSUS
    }*/
    {"filter", filter, &opt.filterWorkflow, "filter chimeric reads",
      "identify reads containing spurious nucleotide order",
      "Annika Jochheim <annika.jochheim@mpibpc.mpg.de>",
      " --reads <fastaFile> [--counts <count.h5>] [--outprefix <string>] [options]",
      FILTER
    },
    {"abundance", abundanceEstimator, &opt.abundanceEstimatorWorkflow, "estimate abundance values",
      "Give for every read an estimated value for the abundance",
      "Annika Jochheim <annika.jochheim@mpibpc.mpg.de>",
      " --reads <fastaFile> [--counts <count.h5>] [--outprefix <string>] [options]",
      ABUNDANCE_ESTIMATOR
    },

    /*** developer tools ***/

    {"profile", profile, &opt.profileWorkflow, "print spaced k-mer count profiles (devtool)",
      "dev tool to write for every sequence the spaced k-mer count profile in a tab separated plain text file",
      "Annika Jochheim <annika.jochheim@mpibpc.mpg.de>",
      " --reads <fastaFile> [--counts <count.h5>] [--outprefix <string>] [options]",
      PROFILE
    },
    {"counts2flat", counts2flat, &opt.counts2flatWorkflow, "print spaced k-mer lookup table (devtool)",
     "transform continous kmers to spaced kmers and print spaced k-mers and counts as stored in lookuptable",
     "Annika Jochheim <annika.jochheim@mpibpc.mpg.de>",
     " --counts <count.h5>] [--outprefix <string>] [options]",
      COUNTS2FLAT
    }
  };

struct Command *getCommand(const char *name) {
  for (size_t idx = 0; idx < commands.size(); idx++) {
    struct Command *command = &commands[idx];
    if (!strcmp(name, command->cmd))
      return command;
  }
  return NULL;
}

void printUsage(const int mode = SIMPLE) {
  std::stringstream usage;
  usage << tool_introduction << "\n\n";
  usage << tool_name << " Version: " << version << "\n";
  usage << "© " << main_author << "\n\n";

  usage << "available commands:\n";
  for (size_t j = 0; j < commands.size(); j++) {
    struct Command &t = commands[j];
    usage << "  " + std::string(t.cmd) << "\t"<< t.descriptShort << "\n";
  }

  //TODO: add extra category for dev tools

  std::cerr << usage.str() << "\n";
}

int main(int argc, const char *argv[]) {

  if (argc < 2) {
    printUsage(SIMPLE);
    return (EXIT_FAILURE);
  }

  if ((argv[1][0] == '-' && argv[1][1] == 'h') ||
      strcmp(argv[1], "--help") == 0) {
    printUsage(EXTENDED);
    return (EXIT_SUCCESS);
  }

  struct Command *command = getCommand(argv[1]);
  if (command == NULL) {
    Info(Info::ERROR) << "Invalid Command: " << argv[1] << "\n";
    printUsage(SIMPLE);
    return (EXIT_FAILURE);
  }

  Info(Info::INFO) <<"\n" << tool_name << " Version: " << version << "\n";
  Info(Info::INFO) << "Execute " << tool_name << " command: " << command->cmd << "\n\n";
  EXIT(command->callerFunction(argc - 2, argv + 2, command));
}
