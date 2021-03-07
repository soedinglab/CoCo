// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include "preprocessing.h"
#include "Lookuptable.h"
#include "HashTable.h"
#include "filehandling.h"
#include "kseq.h"
#include "Info.h"

KSEQ_INIT(int, read)


template<unsigned int LOGINDEX, unsigned int LOGOFFSET>
LookupTableBase *buildLookuptableIter(string countFile, int countMode,
                                      const KmerTranslator &translator,
                                      uint32_t minCount){

  /* get dsk kmer-count storage */
  Storage *storage = StorageFactory(STORAGE_HDF5).load(countFile);
  LOCAL (storage);
  Group &dskGroup = storage->getGroup("dsk");
  string kmerSizeStr = dskGroup.getProperty("kmer_size");
  unsigned int kmerSize = atoi(kmerSizeStr.c_str());

  unsigned int kmerSpan = translator.getSpan();

  if (kmerSize != kmerSpan) {
    Info(Info::ERROR) << "ERROR: kmerSize " << kmerSize << " used in hdf5 file " << countFile.c_str() << " is not supported.\n"
                                                                                                         "Please pre-compute kmer counts with k=41\n";
    return NULL;
  }

  //retrieve the partition holding the couples [kmer,abundance]

  Partition<Count> &solidKmers = dskGroup.getPartition<Count>("solid");
  assert(translator.getSpan() == kmerSize);

  Lookuptable<LOGINDEX,LOGOFFSET> *lookuptable = new Lookuptable<LOGINDEX, LOGOFFSET>(solidKmers.getNbItems(), countMode);
  // fill lookuptable
  {
    Iterator<Count> *it = solidKmers.iterator();
    LOCAL (it);
    typename Kmer<>::ModelCanonical model(kmerSize);

    Info(Info::INFO) << "...construct grids\n";

    // construct grids: count kmers per grid, increase grid values in doing so
    for (it->first(); !it->isDone(); it->next()) {
      const Count &count = it->item();

      packedKmerType packedKmer =
        translator.kmer2minPackedKmer(largeInt2uint128(count.value));
      lookuptable->assignKmertoGrid(packedKmer);
    }

    // shift grid value to grid start positions
    lookuptable->setupIndexGridTable();

    Info(Info::INFO) << "...fill lookuptable\n";
    // add elements, increase grid values in doing so
    for (it->first(); !it->isDone(); it->next()) {
      const Count &count = it->item();
      packedKmerType packedKmer =
        translator.kmer2minPackedKmer(largeInt2uint128(count.value));
      lookuptable->addElement(packedKmer, count.abundance);
    }

    Info(Info::INFO) << "...final setup\n";

    // final setup - shift grid value to grid start positions back
    lookuptable->finalSetupTables(minCount);

    //TODO: sort elements per grid?

    Info(Info::INFO) << "...completed\n";

#ifdef DEBUG
    if(!isValid(*lookuptable, solidKmers, model, translator, minCount))
      return NULL;
#endif
  }
  return lookuptable;
}

LookupTableBase *buildLookuptable(string countFile, int countMode,
                                  const KmerTranslator &translator,
                                  uint32_t minCount) {




  Info(Info::INFO) << "create lookuptable...\n";
  unsigned int logIndexSize, logOffsetSize;
  translator.getBestSplit(logIndexSize, logOffsetSize );
  // create lookuptable
  if(logIndexSize==22 && logOffsetSize==2)
    return buildLookuptableIter<22,2>(countFile, countMode,translator, minCount);
  else if(logIndexSize==24 && logOffsetSize==2)
    return buildLookuptableIter<24,2>(countFile, countMode,translator, minCount);
  else if(logIndexSize==26 && logOffsetSize==2)
    return buildLookuptableIter<26,2>(countFile, countMode,translator, minCount);
  else if(logIndexSize==28 && logOffsetSize==2)
    return buildLookuptableIter<28,2>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==2)
    return buildLookuptableIter<30,2>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==4)
    return buildLookuptableIter<30,4>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==6)
    return buildLookuptableIter<30,6>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==8)
    return buildLookuptableIter<30,8>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==10)
    return buildLookuptableIter<30,10>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==12)
    return buildLookuptableIter<30,12>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==14)
    return buildLookuptableIter<30,14>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==16)
    return buildLookuptableIter<30,16>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==18)
    return buildLookuptableIter<30,18>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==20)
    return buildLookuptableIter<30,20>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==22)
    return buildLookuptableIter<30,22>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==24)
    return buildLookuptableIter<30,24>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==26)
    return buildLookuptableIter<30,26>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==28)
    return buildLookuptableIter<30,28>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==30)
    return buildLookuptableIter<30,30>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==32)
    return buildLookuptableIter<30,32>(countFile, countMode,translator, minCount);
  else if(logIndexSize==30 && logOffsetSize==34)
    return buildLookuptableIter<30,34>(countFile, countMode,translator, minCount);
  else {
    Info(Info::ERROR) << "ERROR: this weight is not implemented, only values 12-32 are supported\n";
    EXIT(EXIT_FAILURE);
  }

}

LookupTableBase *buildHashTable(string seqFile, const KmerTranslator &translator) {

  Info(Info::WARNING) << "WARNING: counting kmers only on seqfile argument. "\
                         "Make sure the reads are not (pre)clustered! \n";

  Info(Info::WARNING) << "WARNING: using internal hash table to count k-mers is not "\
                         "recommended for larger datasets: use --counts instead\n";

  
  unsigned int kmerSpan = translator.getSpan();
  HashTable *hashtable = new HashTable();

  FILE *kmerCountFile = openFileOrDie(seqFile, "r");
  int fd = fileno(kmerCountFile);
  kseq_t *seq = kseq_init(fd);
  spacedKmerType spacedKmer = 0, mask = ((((spacedKmerType) 1) << (spacedKmerType) (kmerSpan * 2)) - 1);
  packedKmerType x;

  Info(Info::INFO) << "count spaced k-mers...\n";
  while (kseq_read(seq) >= 0) {
    const size_t len = seq->seq.l;
    const char *seqNuc = seq->seq.s;
    SeqType seqStr;
    seqStr.reserve(len);

    if (len < kmerSpan)
      continue;

    /* sequence to 2bit representation */
    unsigned int l;
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

  Info(Info::INFO) << "...completed\n";
  return hashtable;
}


//DEBUG
bool isValid(const LookupTableBase &lookuptable,
             Partition<Count> &solidKmers,
             Kmer<>::ModelCanonical &model,
             const KmerTranslator &translator,
             unsigned short threshold) {

  Iterator<Count> *it = solidKmers.iterator();
  LOCAL (it);

  for (it->first(); !it->isDone(); it->next()) {
    const Count &count = it->item();
    packedKmerType packedKmer =
      translator.kmer2minPackedKmer(largeInt2uint128(count.value));

    unsigned int lookupCount = lookuptable.getCount(packedKmer);
    std::string kmer = model.toString(count.value);

    //std::cout << lookupCount << "\t" << count.abundance << std::endl;

    // abundance < threshold:
    // kmer should not exist in lookuptable, except different kmers match to
    // the same packedKmer by using spacemask
    if (count.abundance < threshold && lookupCount != 0 && lookupCount < threshold) {

      Info(Info::ERROR)  << "ERROR: count value " << lookupCount
                         << "found for packed kmer "
                         << packedKmer << "generated from " << kmer
                         << " is smaller than threshold" << threshold << "\n";
      return false;
    }

    // check count values
    if ((unsigned) count.abundance > threshold && (unsigned) count.abundance > lookupCount) {
      Info(Info::ERROR) << "ERROR: check lookuptable counts failed, value " << lookupCount
                        << " for packed kmer " << packedKmer << "generated from"
                        << kmer << "is smaller than count.abundance="
                        << count.abundance << "found in h5 file " << "\n";
      return false;
    }
  }
  return true;
}
