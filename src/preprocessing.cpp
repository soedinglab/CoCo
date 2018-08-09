// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include "preprocessing.h"
#include <gatb/gatb_core.hpp>
#include "Lookuptable.h"
#include "util.h"


bool isValid(const Lookuptable &lookuptable,
              Partition<Count> &solidKmers,
              Kmer<>::ModelCanonical &model,
              const KmerTranslator &translator,
              unsigned short threeshold)
{

  Iterator<Count>* it = solidKmers.iterator();  LOCAL (it);

  for (it->first(); !it->isDone(); it->next())
  {
    const Count& count = it->item();
    kmerType packedKmer =
        translator.kmer2minPackedKmer(largeInt2uint128(count.value));
    //TODO: packedKmer as string
    unsigned int lookupCount = lookuptable.getCount(packedKmer);
    std::string kmer = model.toString(count.value);

    //std::cout << lookupCount << "\t" << count.abundance << std::endl;

    // abundance < threeshold:
    // kmer should not exist in lookuptable, except different kmers match to
    // the same packedKmer by using spacemask
    if (count.abundance < threeshold && lookupCount!=0 && lookupCount < threeshold)
    {

      std::cerr << "ERROR: count value " << lookupCount
                << "found for packed kmer "
                << packedKmer << "generated from " << kmer
                << " is smaller than threeshold" << threeshold << std::endl;
      return false;
    }

    // check count values
    if ((unsigned)count.abundance > threeshold && (unsigned)count.abundance > lookupCount)
    {
      std::cerr << "ERROR: check lookuptable counts failed, value " << lookupCount
                << " for packed kmer " << packedKmer << "generated from"
                << kmer << "is smaller than count.abundance="
                << count.abundance << "found in h5 file " << std::endl;
      return false;
    }
  }
  return true;
}

Lookuptable* buildLookuptable(Storage &storage,
                              const KmerTranslator &translator,
                              size_t minCount,
                              float corrFactor)
{
  //retrieve the partition holding the couples [kmer,abundance]
  Group& dskGroup = storage.getGroup("dsk");
  Partition<Count>& solidKmers = dskGroup.getPartition<Count> ("solid");
  string kmerSizeStr = dskGroup.getProperty ("kmer_size");
  unsigned int kmerSize = atoi(kmerSizeStr.c_str());


  assert(translator.getSpan() == kmerSize);
  fprintf(stderr, "start create lookuptable\n");
  // create lookuptable
  Lookuptable *lookuptable = new Lookuptable(solidKmers.getNbItems(), corrFactor);
  // fill lookuptable
  {
    Iterator<Count>* it = solidKmers.iterator();  LOCAL (it);
    typename Kmer<>::ModelCanonical model(kmerSize);

     fprintf(stderr, "start count kmers - construct grids\n");
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

    fprintf(stderr, "start fill lookuptable\n");
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

    fprintf(stderr, "start final setup for lookuptable\n");
    // shift grid value to grid start positions back
    lookuptable->finalSetupTables(minCount);

    //TODO: sort elements per grid

#ifdef DEBUG
    if(!isValid(*lookuptable,solidKmers,model,translator,minCount))
      return NULL;
#endif
  }
  return lookuptable;
}



