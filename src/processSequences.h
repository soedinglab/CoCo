#ifndef PROCESSREADS
#define PROCESSREADS

#include "CountProfile.h"
#include "Lookuptable.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "util.h"


int processSeqFile(string seqFilename,
                   string resultFilename,
                   const Lookuptable* lookuptable,
                   const KmerTranslator* translator,
                   int (*processCountProfile)(CountProfile &, FILE*),
                   size_t chunkStart=0,
                   size_t chunkEnd=std::numeric_limits<uint64_t>::max());
//TODO: better way for handle chunEnd default?

int processSeqFileParallel(string seqFilename,
                           string resultFilename,
                           const Lookuptable* lookuptable,
                           const KmerTranslator* translator,
                           int (*processCountProfile)(CountProfile &, FILE*),
                           int threadNum);

void process_sampleList(vector<string> *sampleList, std::string resultFileName,
                        Lookuptable *lookuptable,\
                        KmerTranslator *translator,
                        int (*processCountProfile)(CountProfile &, FILE*));

#endif // PROCESSREADS

