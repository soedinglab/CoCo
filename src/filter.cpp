// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <cstdio>
#include <stdio.h>
#include <fcntl.h>
#include <stdexcept>

#include "Command.h"
#include "options.h"
#include "types.h"
#include "KmerTranslator.h"
#include "preprocessing.h"
#include "runner.h"
#include "filehandling.h"

typedef struct {
    FILE *filterReads;
    FILE *cleanedReads;
} FilterArgs;

int filterProcessor(CountProfile &countprofile, void *filterargs)
{
    std::vector<unsigned int> dropPositions = countprofile.getDropPointsInMaximzedProfile();

    bool toFilter = countprofile.checkForRiseAndDropPoints(dropPositions);

    SequenceInfo *seqinfo = countprofile.getSeqInfo();
    if(toFilter)
    {
        fwrite(&seqinfo->sep, sizeof(char), 1, ((FilterArgs *) filterargs)->filterReads);
        fwrite(seqinfo->name.c_str(), sizeof(char), seqinfo->name.size(), ((FilterArgs *) filterargs)->filterReads);
        fwrite(seqinfo->comment.c_str(), sizeof(char), seqinfo->comment.size(), ((FilterArgs *) filterargs)->filterReads);
        fwrite("\n", sizeof(char), 1, ((FilterArgs *) filterargs)->filterReads);
        fwrite(seqinfo->seq.c_str(), sizeof(char), seqinfo->seq.size(), ((FilterArgs *) filterargs)->filterReads);
        fwrite("\n", sizeof(char), 1, ((FilterArgs *) filterargs)->filterReads);
        if (seqinfo->sep == '@')
        {
            fwrite("*", sizeof(char), 1, ((FilterArgs *) filterargs)->filterReads);
            fwrite(seqinfo->qual.c_str(), sizeof(char), seqinfo->qual.size(), ((FilterArgs *) filterargs)->filterReads);
            fwrite("\n", sizeof(char), 1, ((FilterArgs *) filterargs)->filterReads);
        }
    }
    else
    {
        fwrite(&seqinfo->sep, sizeof(char), 1, ((FilterArgs *) filterargs)->cleanedReads);
        fwrite(seqinfo->name.c_str(), sizeof(char), seqinfo->name.size(), ((FilterArgs *) filterargs)->cleanedReads);
        fwrite(seqinfo->comment.c_str(), sizeof(char), seqinfo->comment.size(), ((FilterArgs *) filterargs)->cleanedReads);
        fwrite("\n", sizeof(char), 1, ((FilterArgs *) filterargs)->cleanedReads);
        fwrite(seqinfo->seq.c_str(), sizeof(char), seqinfo->seq.size(), ((FilterArgs *) filterargs)->cleanedReads);
        fwrite("\n", sizeof(char), 1, ((FilterArgs *) filterargs)->cleanedReads);
        if (seqinfo->sep == '@')
        {
            fwrite("*", sizeof(char), 1, ((FilterArgs *) filterargs)->cleanedReads);
            fwrite(seqinfo->qual.c_str(), sizeof(char), seqinfo->qual.size(), ((FilterArgs *) filterargs)->cleanedReads);
            fwrite("\n", sizeof(char), 1, ((FilterArgs *) filterargs)->cleanedReads);
        }
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

    string outprefix;
    if(opt.OP_OUTPREFIX.isSet)
        outprefix = opt.outprefix;
    else
        outprefix = getFilename(seqFile);
    string ext = getFileExtension(seqFile);
    FilterArgs filterargs = {openFileOrDie(outprefix + ".filtered" + ext, "w"),
                             openFileOrDie(outprefix + ".cleaned" + ext, "w")};

    processSeqFile(seqFile, lookuptable, translator, filterProcessor, &filterargs);

    fclose(filterargs.filterReads);
    fclose(filterargs.cleanedReads);

    delete lookuptable;
    delete translator;

    return EXIT_SUCCESS;
}


