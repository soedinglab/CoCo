// Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>

#ifndef RUNNER
#define RUNNER

#include "CountProfile.h"
#include "LookuptableBase.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "util.h"


int processReads(string readsname,
                 LookupTableBase *lookuptable,
                 const KmerTranslator *translator,
                 int (*processCountProfile)(CountProfile &, void *, bool),
                 void *processArgs,
                 int skip,
                 bool silent = false,
                 size_t chunkStart = 0,
                 size_t chunkEnd = std::numeric_limits<uint64_t>::max());

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

