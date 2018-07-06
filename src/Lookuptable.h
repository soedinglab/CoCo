// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de> and
// Martin Steinegger <martin.steinegger@mpibpc.mpg.de>

#ifndef LOOKUPTABLE_H
#define LOOKUPTABLE_H

#include <cstdlib>
#include "mathsupport.h"
#include "kmer.h"

/* k-mer-index = (first p bits << I) | (last I bits) */
#define LOGINDEXSIZE 30 /* p bits, gurantee ~8GB for indexGridTable*/
#define LOGOFFSETSIZE 24 /* I bits */

class Lookuptable
{

private:

  size_t *indexGridTable;
  size_t indexGridTableSize = ipow(2,LOGINDEXSIZE)+1;
  struct __attribute__((__packed__)) IndexEntry {
    unsigned int indexOffset : LOGOFFSETSIZE;
    unsigned int count;
    static bool compareByKmer(IndexEntry first, IndexEntry second){
      return (first.indexOffset < second.indexOffset );
    }
  };
  IndexEntry * offsetTable;
  size_t numberItems;
  size_t maxNumberItems;
  //TODO: size_t countThreeshold

  size_t _offsetmask;
  size_t _indexmask;

  inline size_t getGridPosition(kmerType kmer) const;
  inline size_t getOffset(kmerType kmer) const;

public:

  Lookuptable(){};

  Lookuptable(const size_t nbItems);

  ~Lookuptable();

  void assignKmertoGrid(kmerType kmer);
  void setupIndexGridTable();
  size_t addElement(kmerType kmer, unsigned int count);
  void finalSetupTables(size_t countThreeshold=1);
};

#endif // LOOKUPTABLE_H
