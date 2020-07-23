#include "preprocessing.h"
#include "Lookuptable.h"
#include "HashTable.h"
#include "filehandling.h"
#include "kseq.h"
KSEQ_INIT(int, read)


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

LookupTableBase* buildLookuptable(string countFile,
                              const KmerTranslator &translator,
                              size_t minCount,
                              float corrFactor)
{
    unsigned int kmerSpan = translator.getSpan();

    /* get dsk kmer-count storage */
    Storage* storage = StorageFactory(STORAGE_HDF5).load(countFile);
    LOCAL (storage);
    Group& dskGroup = storage->getGroup("dsk");
    string kmerSizeStr = dskGroup.getProperty ("kmer_size");
    unsigned int kmerSize = atoi(kmerSizeStr.c_str());

    if (kmerSize != kmerSpan)
    {
        fprintf(stderr, "kmerSize %u used in hdf5 file %s is not supported.\n"
                        "Please precompute kmer counts with k=41\n", kmerSize, countFile.c_str());
        return NULL;
    }

  //retrieve the partition holding the couples [kmer,abundance]

  Partition<Count>& solidKmers = dskGroup.getPartition<Count> ("solid");

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

LookupTableBase *buildHashTable(string seqFile, const KmerTranslator &translator)
{
    unsigned int kmerSpan = translator.getSpan();
    HashTable *hashtable = new HashTable();

    FILE *kmerCountFile = openFileOrDie(seqFile, "r");
    int fd = fileno(kmerCountFile);
    kseq_t *seq = kseq_init(fd);
    spacedKmerType spacedKmer, mask = ((((spacedKmerType) 1) << (spacedKmerType)( kmerSpan * 2)) - 1);
    kmerType x;
    while (kseq_read(seq) >= 0) {
        const size_t len = seq->seq.l;
        const char *seqNuc = seq->seq.s;
        const char *seqName = seq->name.s;
        SeqType seqStr;
        seqStr.reserve(len);

        if (len < kmerSpan)
            continue;

        /* sequence to 2bit representation */
        int l;
        for (unsigned int pos = l = 0; pos < len; pos++) {
            int c = res2int[(int) seqNuc[pos]];
            if (c != -1) {
                spacedKmer = (spacedKmer << 2 | c) & mask;
                if (++l >= kmerSpan) {
                    x = translator.kmer2minPackedKmer(spacedKmer);
                    hashtable->increaseCount(x);
                }
            } else {
                l = 0;
                spacedKmer = 0;
                x = 0;
            }
        }
        //TODO: add kmers with non informative N's ?
    }
    kseq_destroy(seq);
    fclose(kmerCountFile);
    return hashtable;
}



