#ifndef PREPROCESSING
#define PREPROCESSING

// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <gatb/gatb_core.hpp>
#include "Lookuptable.h"

Lookuptable* buildLookuptable(const Storage &storage,
                              const KmerTranslator & translator, size_t minCount)
{

  LOCAL (storage);
  string kmerSizeStr = storage->getGroup("dsk").getProperty ("kmer_size");
  unsigned long kmerSize = atol(kmerSizeStr.c_str());

  //retrieve the partition holding the couples [kmer,abundance]
  Group& dskGroup = storage->getGroup("dsk");
  Partition<Count>& solidKmers = dskGroup.getPartition<Count> ("solid");

  unsigned short span = translator.getSpan();
  unsigned short weight = translator.getWeight();

  assert(span == kmerSize);

  // create lookuptable
  Lookuptable *lookuptable = new Lookuptable(kmerSize, solidKmers.getNbItems());
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

    //TODO: resize always necessary? trade of memory and time?
    // shift grid value to grid start positions back
    lookuptable->finalSetupTables(true, minCount);

    //TODO: sort elements per grid

    //TODO: valid check
  }
  return lookuptable;
}

#endif // PREPROCESSING

