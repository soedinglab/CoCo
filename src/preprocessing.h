#ifndef PREPROCESSING
#define PREPROCESSING

// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <gatb/gatb_core.hpp>
#include "Lookuptable.h"

typedef Kmer<>::Count Count;

Lookuptable* buildLookuptable(Storage &storage,
                              const KmerTranslator &translator,
                              size_t minCount)
{
  //retrieve the partition holding the couples [kmer,abundance]
  Group& dskGroup = storage.getGroup("dsk");
  Partition<Count>& solidKmers = dskGroup.getPartition<Count> ("solid");
  string kmerSizeStr = dskGroup.getProperty ("kmer_size");
  unsigned int kmerSize = atoi(kmerSizeStr.c_str());

  unsigned short span = translator.getSpan();
  unsigned short weight = translator.getWeight();

  assert(span == kmerSize);

  // create lookuptable
  Lookuptable *lookuptable = new Lookuptable(solidKmers.getNbItems());
  // fill lookuptable
  {
    Iterator<Count>* it = solidKmers.iterator();  LOCAL (it);
    typename Kmer<>::ModelCanonical model(kmerSize);

    // construct grids: count kmers per grid, increase grid values in doing so
    for (it->first(); !it->isDone(); it->next())
    {
      const Count& count = it->item();

      kmerType packedKmer =
          translator.kmer2minPackedKmer(largeInt2uint128(count.value));
      lookuptable->assignKmertoGrid(packedKmer);
    }

    // shift grid value to grid start positions
    lookuptable->setupIndexGridTable();

    // add elements, increase grid values in doing so
    for (it->first(); !it->isDone(); it->next())
    {
      const Count& count = it->item();
      //std::string kmer = model.toString(count.value);

      kmerType packedKmer =
          translator.kmer2minPackedKmer(largeInt2uint128(count.value));
      //TODO: duplicate code, better solution?

      lookuptable->addElement(packedKmer, count.abundance);
    }

    // shift grid value to grid start positions back
    lookuptable->finalSetupTables(minCount);

    //TODO: sort elements per grid

    //TODO: valid check
  }
  return lookuptable;
}

#endif // PREPROCESSING

