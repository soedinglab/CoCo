// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de> and
// Martin Steinegger <martin.steinegger@mpibpc.mpg.de>

#ifndef LOOKUPTABLE_H
#define LOOKUPTABLE_H

#include <cstdlib>
#include "mathsupport.h"

/* k-mer-index = (first p bits << I) | (last I bits) */
#define LOGINDEXSIZE 30 /* p bits, gurantee 8GB for indexGridTable*/
#define LOGOFFSETSIZE 24 /* I bits */

class Lookuptable
{

private:

  size_t *indexGridTable;
  size_t indexGridTableSize = ipow(2,LOGINDEXSIZE);
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

public:

  Lookuptable(){};

  Lookuptable(const size_t kmerWeight,
              const size_t nbItems);
};

#endif // LOOKUPTABLE_H
