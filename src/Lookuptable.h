// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de> and
// Martin Steinegger <themartinsteinegger@gmail.com>

#ifndef LOOKUPTABLE_H
#define LOOKUPTABLE_H

#include <cstdlib>
#include "mathsupport.h"
#include "LookuptableBase.h"
#include "kmer.h"


/* k-mer-index = (first p bits << I) | (last I bits) */
//#define LOGINDEXSIZE 30 /* p bits, guarantee ~8GB for indexGridTable*/
//#define LOGOFFSETSIZE 24 /* I bits */

template<unsigned int LOGINDEXSIZE, unsigned int LOGOFFSETSIZE>
class Lookuptable : public LookupTableBase {

private:

  size_t *indexGridTable;
  size_t indexGridTableSize = ipow(2, LOGINDEXSIZE) + 1;

  struct __attribute__((__packed__)) IndexEntry {
    unsigned long long indexOffset : LOGOFFSETSIZE;
    unsigned int count;
/*
    static bool compareByKmer(IndexEntry first, IndexEntry second) {
      return (first.indexOffset < second.indexOffset);
    }*/
  };

  IndexEntry *offsetTable;
  size_t numberItems;
  size_t maxNumberItems;

  size_t _offsetmask;
  size_t _indexmask;

  int mode;

  inline size_t getGridPosition(packedKmerType kmer) const;

  inline size_t getOffset(packedKmerType kmer) const;

  inline std::pair<size_t, size_t> getIndexGridRange(packedKmerType kmer) const;

public:

  Lookuptable() = default;

  Lookuptable(size_t nbItems, int countMode);

  ~Lookuptable();

  void assignKmertoGrid(packedKmerType kmer);

  void setupIndexGridTable();

  size_t addElement(packedKmerType kmer, unsigned int count);

  void finalSetupTables(size_t countThreshold = 1);

  unsigned int getCount(const packedKmerType kmer) const;

  void iterateOverAll(FILE* fp) const;

  bool decreaseCount(packedKmerType kmer);

  bool increaseCount(packedKmerType kmer);
  void printSize(){
    std::cout << sizeof(IndexEntry) << std::endl;
    std::cout << "IndexTable: " << sizeof(*indexGridTable)*indexGridTableSize << std::endl;
    std::cout << "OffsetTable: " << sizeof(*offsetTable)*maxNumberItems << std::endl;
    //exit(1);
  }

};

#endif // LOOKUPTABLE_H
