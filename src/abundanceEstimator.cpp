// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include <cstdio>
#include <fcntl.h>
#include <stdexcept>
#include <thread>
#include <vector>

#include "Command.h"
#include "options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "Lookuptable.h"
#include "CountProfile.h"
#include "preprocessing.h"
#include "runner.h"
#include "filehandling.h"

typedef struct {
    FILE *abundanceFile;
} AbundanceEstimatorArgs;


int abundanceEstimatationProcessor(CountProfile &countprofile, void *abundanceargs)
{
  // estimate abundance value
  // TODO
}


int abundanceEstimator(int argc, const char **argv, const Command* tool)
{
    int exit_code=0;
    Options &opt = Options::getInstance();
    opt.parseOptions(argc, argv, *tool);

    //TODO: print parameters

    // TODO:check parameter and if files exists

    initialize();
    KmerTranslator *translator = new KmerTranslator();
    string seqFile = opt.seqFile;


    LookupTableBase *lookuptable;

    // use precomputed counts and fill lookuptable
    if (opt.OP_COUNT_FILE.isSet) {

        string countFile = opt.countFile;
        lookuptable = buildLookuptable(countFile, *translator, 0, 1);
    }
    else { // count k-mers itself and fill hash-lookuptable

        lookuptable = buildHashTable(seqFile, *translator);
    }

    if (lookuptable == NULL) {

        fprintf(stderr,"Generating lookuptablefailed\n");
        return EXIT_FAILURE;
    }

    AbundanceEstimatorArgs abundanceargs = {openFileOrDie("abundance", "w")};

    if (opt.threads == 1)
    {
      exit_code = processSeqFile(seqFile, lookuptable, translator, abundanceEstimatationProcessor, &abundanceargs);
      if (exit_code != 0)
      {
        std::cerr << "ERROR processing sequence file " << seqFile << std::endl;
      }
    }

    //    else
//    {
//      retval = processSeqFileParallel(seqFile,
//                                      resultFile,
//                                      lookuptable,
//                                      translator,
//                                      writeAbundanceEstimation,
//                                      opt.threads);
//    }
//  }
//

  fclose(abundanceargs.abundanceFile);
  delete lookuptable;
  delete translator;

  return exit_code;
}
