/*
 * Copyright (C) 2018 Annika Seidel <annika.seidel@mpibpc.mpg.de>
 *
 */

#include <cstdio>
#include <fcntl.h>
#include <stdexcept>
#include <vector>
#include <gatb/gatb_core.hpp>

#include "Command.h"
#include "options.h"
#include "types.h"
#include "KmerTranslator.h"

#include "CountProfile.h"
#include "preprocessing.h"
#include "processSequences.h"
#include "filehandling.h"


/*int showProfile(CountProfile &countprofile, FILE *resultFile)
{
  fprintf(resultFile,"#%s\n", countprofile.getSeqName());
  countprofile.showProfile(resultFile);

}*/

int profile(int argc, const char **argv, const Command* tool) {
//    /*Options &opt = Options::getInstance();
//    opt.parseOptions(argc, argv, *tool);
//    printf("seqFile: %s\n", opt.seqFile.c_str());
//    printf("kmerCounFile: %s\n", opt.kcFile.c_str());
//    initialize();
//    int span = 41;
//    int weight = 27;
//    const char *seqFilename = opt.seqFile.c_str();
//    const char * kmerCountFileName = opt.kcFile.c_str();
//    KmerTranslator *translator = new KmerTranslator(span, weight);
//    HashTable *lookuptable = new HashTable();
//    fprintf(stderr, "Counting kmers...\n");
//    {// Fill hashlookuptable
//        FILE *kmerCountFile = openFileOrDie(kmerCountFileName, "r");
//        int fd = fileno(kmerCountFile);
//        kseq_t *seq = kseq_init(fd);
//        spacedKmerType spacedKmer, mask = ((((spacedKmerType) 1) << (spacedKmerType)( span * 2)) - 1);
//        kmerType x;
//        while (kseq_read(seq) >= 0) {
//            const size_t len = seq->seq.l;
//            const char *seqNuc = seq->seq.s;
//            const char *seqName = seq->name.s;
//            SeqType seqStr;
//            seqStr.reserve(len);
//
//            unsigned int kmerSpan = translator->getSpan();
//            if (len < kmerSpan) {
//                fprintf(stderr, "WARNING: sequence %s is too short, it'll be skipped\n", seqName);
//                continue;
//            }
//
//            *//* sequence to 2bit representation *//*
//            int l;
//            for (unsigned int pos = l = 0; pos < len; pos++) {
//                int c = res2int[(int) seqNuc[pos]];
//                if (c != -1) {
//                    spacedKmer = (spacedKmer << 2 | c) & mask;                  // forward strand
//                    if (++l >= span) {
//                        x = translator->kmer2minPackedKmer(spacedKmer);
//                        lookuptable->increaseCount(x);
//                    }
//                } else {
//                    l = 0;
//                    spacedKmer = 0;
//                    x = 0;
//                }
//            }
//            //TODO: add kmers with non informative N's ?
//        }
//        kseq_destroy(seq);
//        fclose(kmerCountFile);
//    }
//    fprintf(stderr, "...done\nProcessing Sequences...");
//    //processSeqFile(seqFilename, "resultfile.txt", lookuptable, translator,showProfile);*/
    return EXIT_SUCCESS;
}

