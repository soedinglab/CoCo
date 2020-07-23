#include <cstdio>
#include <fcntl.h>
#include <stdexcept>

#include "Command.h"
#include "options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "preprocessing.h"
#include "processSequences.h"
#include "filehandling.h"

typedef struct {
    FILE *filterReads;
    FILE *keepReads;
} FilterArgs;

int filterProcessor(CountProfile &countprofile, void *filterargs)
{

    std::vector<unsigned int> dropPositions = countprofile.getDropPointsInMaximzedProfile();

    bool toFilter = countprofile.checkForRiseAndDropPoints(dropPositions);
    if(toFilter)
    {
        fprintf(((FilterArgs *)filterargs)->filterReads, "%s\n", countprofile.getSeqName());
    }

}

int filter(int argc, const char **argv, const Command* tool)
{
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

    FilterArgs filterargs = {openFileOrDie("filtered", "w"),
                             openFileOrDie("keep", "w")};

    processSeqFile(seqFile, lookuptable, translator, filterProcessor, &filterargs);

    fclose(filterargs.filterReads);
    fclose(filterargs.keepReads);

    /*
     *string filename = kmerCountFile;
    size_t lastdot = filename.find_last_of(".");
    if (lastdot != std::string::npos)
        filename=filename.substr(0, lastdot);
    string resultFilename=(string("count_profile.")+string(basename(filename.c_str()))+string(".txt"));*/

    delete lookuptable;
    delete translator;

    return EXIT_SUCCESS;
}


