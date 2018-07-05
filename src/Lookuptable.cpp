#include <assert.h>
#include <string.h>
#include "Lookuptable.h"

Lookuptable::Lookuptable(const size_t kmerWeight, const size_t nbItems)
{
  assert(LOGINDEXSIZE+LOGOFFSETSIZE == 2*kmerWeight);

  //TODO: store one more grid entry?
  indexGridTable = (size_t *) calloc(indexGridTableSize, sizeof(size_t));
  offsetTable = (IndexEntry *) calloc(nbItems, sizeof(IndexEntry));
  maxNumberItems = nbItems;
  numberItems = 0;

  //_indexmask = (ipow(2,LOGINDEXSIZE)-1) << (2*kmerSize-LOGINDEXSIZE);
  //_offsetmask = ipow(2,LOGOFFSETSIZE)-1
}

