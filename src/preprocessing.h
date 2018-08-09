#ifndef PREPROCESSING
#define PREPROCESSING

// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <gatb/gatb_core.hpp>
#include "Lookuptable.h"
#include "KmerTranslator.h"
#include "util.h"

typedef Kmer<>::Count Count;

bool isValid(const Lookuptable &lookuptable,
              Partition<Count> &solidKmers,
              Kmer<>::ModelCanonical &model,
              const KmerTranslator &translator,
              unsigned short threeshold);

Lookuptable* buildLookuptable(Storage &storage,
                              const KmerTranslator &translator,
                              size_t minCount,
                              float corrFactor);


#endif // PREPROCESSING

