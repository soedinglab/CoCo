// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#ifndef PREPROCESSING
#define PREPROCESSING


#include <gatb/gatb_core.hpp>
#include "LookuptableBase.h"
#include "KmerTranslator.h"
#include "util.h"

typedef Kmer<>::Count Count;


LookupTableBase *buildLookuptable(string countFile, int countMode,
                                  const KmerTranslator &translator,
                                  uint32_t minCount);

LookupTableBase *buildHashTable(string seqFile, const KmerTranslator &translator);

/*bool isValid(const Lookuptable &lookuptable,
              Partition<Count> &solidKmers,
              Kmer<>::ModelCanonical &model,
              const KmerTranslator &translator,
              unsigned short threeshold);
*/

#endif // PREPROCESSING

