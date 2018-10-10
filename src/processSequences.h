#ifndef PROCESSREADS
#define PROCESSREADS

#include "CountProfile.h"
#include "Lookuptable.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "util.h"


int processSeqFile(const char* seqFilename,
                   const char* resultFilename,
                   const Lookuptable &lookuptable,
                   const KmerTranslator &translator,
                   int (*processCountProfile)(CountProfile &, FILE*));

int processSeqFile(string seqFilename,
                   string resultFilename,
                   const Lookuptable &lookuptable,
                   const KmerTranslator &translator,
                   int (*processCountProfile)(CountProfile &, FILE*));

int processSeqFileParallel(const char* seqFilename, int threadNum);

void process_sampleList(vector<string> *sampleList, std::string resultFileName,
                        Lookuptable *lookuptable,\
                        KmerTranslator *translator,
                        int (*processCountProfile)(CountProfile &, FILE*));

#endif // PROCESSREADS

