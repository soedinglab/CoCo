//Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>


#include "src/util.h"
#include "src/options.h"
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


const char *tool_binary = "coco";
const char *tool_name = "CoCo";
const char *tool_introduction = "CoCo is an open-source software suite for "\
                                "different COnsensus COmputation applications "\
                                "on short reads and contigs";

const char *main_author = "Annika Seidel (annika.seidel@mpibpc.mpg.de)";

extern int profile(int argc, const char **argv, const struct Command *tool);

extern int filter(int argc, const char **argv, const struct Command *tool);

extern int abundanceEstimator(int argc, const char **argv, const struct Command *tool);

extern int consensus(int argc, const char **argv, const struct Command *tool);

Options &opt = Options::getInstance();
std::vector<struct Command> commands =
  {
    {"profile", profile, &opt.profileWorkflow, "write spaced k-mer count profiles (devtool)",
      "dev tool to write for every read (contig) the spaced k-mer count profile in a tab separated plain text file",
      "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
      " --seqFile <fastaFile> [--count <count.h5>] [--outprefix <string>]",
      PROFILE
    },
    {"filter", filter, &opt.filterWorkflow, "denoise read file",
      "identify reads containing remarkably/irregular nucleotide sequence as chimeras, indels, ... ",
      "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
      " --seqFile <fastaFile> [--count <count.h5>] [--outprefix <string>]",
      FILTER
    },
    {"abundance", abundanceEstimator, &opt.abundanceEstimatorWorkflow,
                                                              "estimate abundance values",
      "Gives for every read an estimated value for the abundance",
      "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
      " --seqFile <fastaFile> [--count <count.h5>] [--outprefix <string>]",
      ABUNDANCE_ESTIMATOR
    },

    {"consensus", consensus, &opt.consensusWorkflow, "calculate consensus reads ",
      "calculate for every read the consensus nucleotide sequence",
      "Annika Seidel <annika.seidel@mpibpc.mpg.de>",
      "--seqFile <fastaFile> [--count <count.h5>] [--outprefix <string>]",
      CONSENSUS
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
    usage << "  " + std::string(t.cmd) << "\t"\
 << t.descriptShort << "\n";
  }

  //TODO: add EXTENDED mode for dev tools, show when -h/--help is used

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
    fprintf(stderr, "Invalid Command: %s\n", argv[1]);
    printUsage(SIMPLE);
    return (EXIT_FAILURE);
  }

  Info(Info::ERROR) << "execute " << tool_name << " tool: " << command->cmd << "\n";
  EXIT(command->callerFunction(argc - 1, argv + 1, command));
}
