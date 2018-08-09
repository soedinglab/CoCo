#ifndef PROCESSREADS
#define PROCESSREADS

#include "CountProfile.h"
#include "Lookuptable.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "util.h"



FILE* openFileOrDie(std::string fileName, const char * mode);

FILE* openFileOrDie(const char *fileName, const char * mode);


int processReadFile(const char* readFilename,
                    const char* resultFilename,
                    const Lookuptable &lookuptable,
                    const KmerTranslator &translator,
                    int (*processCountProfile)(CountProfile &, FILE*));

int processReadFile(string readFilename,
                    string resultFilename,
                    const Lookuptable &lookuptable,
                    const KmerTranslator &translator,
                    int (*processCountProfile)(CountProfile &, FILE*));

void process_sampleList(vector<string> *sampleList, std::string resultFileName,
                        Lookuptable *lookuptable,\
                        KmerTranslator *translator,
                        int (*processCountProfile)(CountProfile &, FILE*));

#endif // PROCESSREADS

