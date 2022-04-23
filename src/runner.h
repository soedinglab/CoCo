// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>

#ifndef RUNNER
#define RUNNER

#include "CountProfile.h"
#include "LookuptableBase.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "util.h"


int processReads(const string &readsname,
                 LookupTableBase *lookuptable,
                 const KmerTranslator *translator,
                 int (*processCountProfile)(CountProfile &, void *, bool),
                 void *processArgs,
                 int skip,
                 bool silent = false);

int processPairedReads(const string &forwardReads,
                       const string &reverseReads,
                       LookupTableBase *lookuptable,
                       const KmerTranslator *translator,
                       int (*processCountProfile)(CountProfile &, CountProfile &, void *, bool),
                       void *processArgs,
                       int skip);

/*int processreadsParallel(string readsname,
                           string resultFilename,
                           const LookupTableBase* lookuptable,
                           const KmerTranslator* translator,
                           int (*processCountProfile)(CountProfile &, FILE*),
                           int threadNum);

void process_sampleList(vector<string> *sampleList, std::string resultFileName,
                        Lookuptable *lookuptable,\
                        KmerTranslator *translator,
                        int (*processCountProfile)(CountProfile &, FILE*));*/

#endif // RUNNER

