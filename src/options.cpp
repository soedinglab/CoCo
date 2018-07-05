#include <iostream>
#include <sstream>
#include <string>
#include <getopt.h>
#include "options.h"
#include "util.h"

Options* Options::instance = NULL;

Options::Options():
  OP_SAMPLE_LIST(OP_SAMPLE_LIST_ID,"sampleList", "--sampleList",
  "name of the file, which contains a list of fasta files line by line",
  typeid(std::string), (void *) &sampleListFile)
  ,
  OP_KC_LIST(OP_KC_LIST_ID,"kmerCountList", "--kmerCountList ",
  "name of the file, which contains k-mer/count files line by line, format: hdf5",
  typeid(std::string),  (void *) &kmerCountListFile)
  ,
  OP_READ_AVERAGELEN_LIST(OP_READ_AVERAGELEN_LIST_ID,"readAvgLenList",
  "--readAvgLenList", "name of the file, which contains a list of the "\
  "average readlength (consider only reads with length > k) per k-mer/count "\
  "file", typeid(std::string),  (void *) &readAvgLenFile)
  ,
  OP_KMER_WEIGHT(OP_KMER_WEIGHT_ID,"kmerWeight",
  "--kmerWeight", "number of informative positions in a k-mer pattern, "
  "default: 27", typeid(int),  (void *) &kmerWeight)
{
  if (instance)
  {
    std::cerr << "Parameter instance already exists!\n";
    abort();
  }
  instance = this;

  //pcoverage
  pcoverageWorkflow.push_back(OP_SAMPLE_LIST);
  pcoverageWorkflow.push_back(OP_KC_LIST);
  pcoverageWorkflow.push_back(OP_READ_AVERAGELEN_LIST);
  pcoverageWorkflow.push_back(OP_KMER_WEIGHT);

  setDefaults();
}


void Options::setDefaults()
{
  sampleListFile = NULL;
  kmerCountListFile = NULL;
  readAvgLenFile = NULL;
  kmerWeight = 27;
}

void printToolUsage(const ToolInfo &tool, const int FLAG)
{
  std::stringstream usage;
  usage << "coco " << tool.cmd << "\n";
  usage << tool.descriptLong << "\n\n";
  usage << "© " << tool.author << "\n\n";

  usage << "Usage: coco " << tool.cmd << " " << tool.usage << "\n\n";

  if (FLAG == EXTENDED)
  {
    const std::vector<cocoOption>& options = *tool.opt;
    for (size_t idx = 0; idx < options.size(); idx++)
    {
      //TODO: layout, maxwidth + 2 space in the beginning
      usage << options[idx].display << "\t" << options[idx].description << "\n";
    }
  }
  std::cerr << usage.str() << "\n";
}

void Options::parseOptions(int argc, const char *argv[],
                           const ToolInfo& tool)
{
  std::vector<cocoOption>& options = *tool.opt;
  int opt, longIndex = 0;
  extern char *optarg;

  static const struct option longOpts[] = {
      {"help", no_argument, NULL, 'h' },
      {"sampleList", required_argument, NULL, 0},
      {"kmerCountList", required_argument, NULL, 0},
      {"readAvgLenList", required_argument, NULL, 0},
      {"kmerWeight", required_argument, NULL, 0},
      { NULL, no_argument, NULL, 0 }
  };

  while ((opt = getopt_long(argc, (char**)argv, "h", longOpts, &longIndex )) != -1)
  {
     switch(opt)
     {
        case 'h':
          printToolUsage(tool, EXTENDED);
          EXIT(EXIT_SUCCESS);
        case 0:
        {
          std::string optname(longOpts[longIndex].name);
          for(size_t idx = 0; idx < options.size(); idx++)
          {
             if (optname.compare("--help") == 0)
             {
                printToolUsage(tool, EXTENDED);
                EXIT(EXIT_SUCCESS);
             }

             if(optname.compare(options[idx].name) == 0)
             {
               if (options[idx].isSet)
               {
                 fprintf(stderr, "Duplicate option %s\n\n",options[idx].display);
                 printToolUsage(tool,SIMPLE);
                 EXIT(EXIT_FAILURE);
               }

               if (typeid(std::string) == options[idx].type)
               {
                 if(optarg[0] == '-')
                 {
                   fprintf(stderr, "Invalid argument %s for option %s\n",
                                   optarg, options[idx].display);
                   EXIT(EXIT_FAILURE);
                 }
                 std::string val(optarg);
                 if(val.length() != 0)
                 {
                   std::string * currVal = ((std::string *)options[idx].value);
                   currVal->assign( val );
                   options[idx].isSet = true;
                 }
               }
               else if (typeid(int) == options[idx].type)
               {
                 int val = atoi(optarg);
                 if (val == 0)
                 {
                   fprintf(stderr, "Invalid argument %s for option %s\n",
                                   optarg, options[idx].display);
                   EXIT(EXIT_FAILURE);
                 }
                 *((int *) options[idx].value) = val;
                 options[idx].isSet = true;

               }
               else
               {
                 fprintf(stderr, "Wrong option type in parseOptions. "\
                         "please send an error report to the developers.");
                 EXIT(EXIT_FAILURE);
               }
               continue;
             }
           }
          //assert(false); //TODO
          break;
        }

       case '?':
         /* error message already printed by getopt function*/
         printToolUsage(tool, SIMPLE);
         EXIT(EXIT_FAILURE);
       default:
         abort ();
    }
  }

  //TODO: check parameter count, superflous parameter?
  //TODO: check requiered parameter
  //TODO: flag für print Parameter
}