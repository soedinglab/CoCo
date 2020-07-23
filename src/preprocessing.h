#ifndef PREPROCESSING
#define PREPROCESSING


#include <gatb/gatb_core.hpp>
#include "LookuptableBase.h"
#include "KmerTranslator.h"
#include "util.h"

typedef Kmer<>::Count Count;

/*bool isValid(const Lookuptable &lookuptable,
              Partition<Count> &solidKmers,
              Kmer<>::ModelCanonical &model,
              const KmerTranslator &translator,
              unsigned short threeshold);
*/
LookupTableBase* buildLookuptable(string countFile,
                              const KmerTranslator &translator,
                              size_t minCount,
                              float corrFactor);

LookupTableBase* buildHashTable(string seqFile, const KmerTranslator &translator);


#endif // PREPROCESSING

