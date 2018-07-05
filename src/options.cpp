#include <iostream>
#include <getopt.h>
#include "options.h"
#include "util.h"

Options* Options::instance = NULL;

Options::Options():
  OP_SAMPLE_LIST(OP_SAMPLE_LIST_ID,"sampleList", "sample list",
  "name of the file, which contains a list of fasta files line by line",
  typeid(std::string), (void *) &sampleListFile)
  ,
  OP_KC_LIST(OP_KC_LIST_ID,"kmerCountList", "kmer count list ",
  "name of the file, which contains k-mer/count files line by line, format: hdf5",
  typeid(std::string),  (void *) &kmerCountListFile)
{
  if (instance) {
    std::cerr << "Parameter instance already exists!\n";
    abort();
  }
  instance = this;

  //pcoverage
  pcoverageWorkflow.push_back(OP_SAMPLE_LIST);
  pcoverageWorkflow.push_back(OP_KC_LIST);

  setDefaults();
}


void Options::setDefaults()
{
  sampleListFile = NULL;
  kmerCountListFile = NULL;
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
      { NULL, no_argument, NULL, 0 }
  };

  while ((opt = getopt_long(argc, (char**)argv, "h", longOpts, &longIndex )) != -1)
  {
     switch(opt)
     {
        case 'h':
          //TODO: usage
          printf("usage\n");
          EXIT(EXIT_SUCCESS);
        case 0:
        {
          std::string optname(longOpts[longIndex].name);
          for(size_t idx = 0; idx < options.size(); idx++)
          {
             if (optname.compare("--help") == 0)
             {
                //TODO: usage
                EXIT(EXIT_SUCCESS);
             }

             if(optname.compare(options[idx].name) == 0)
             {
               //TODO: is_set?
               printf("optarg: %s\n", optarg);
               if (typeid(std::string) == options[idx].type)
               {
                 std::string val(optarg);
                 if(val.length() != 0)
                 {
                   std::string * currVal = ((std::string *)options[idx].value);
                   currVal->assign( val );
                   //TODO: set is_set true
                   //TODO: check start with "-"
                 }
               }
               else
               {
                 fprintf(stderr, "Wrong parameter type in parseOptions. "\
                         "please send an error report to the developers.");
                 EXIT(EXIT_FAILURE);
               }
               continue;
               //TODO: check and store optarg value
             }
           }
          //assert(false); //TODO
          break;
        }

       case '?':
         /* error message already printed by getopt function*/
         EXIT(EXIT_FAILURE);
       default:
         abort ();
    }
  }
  //TODO: check parameter count
  //TODO: flag für print Parameter
}
