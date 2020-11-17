// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#ifndef HASHTABLE_H
#define HASHTABLE_H

#include "LookuptableBase.h"
#include "khashl.h" // hash table

KHASHL_MAP_INIT(, kc_c1_t, kc_c1, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)

class HashTable : public LookupTableBase {

private:
  kc_c1_t *hashTable;
public:
  HashTable() {
    hashTable = kc_c1_init();
  }

  unsigned int getCount(const packedKmerType kmer) const {
    khint_t itr;
    int absent;
    itr = kc_c1_put(this->hashTable, kmer, &absent);
    if (absent) return 1;
    return kh_val(this->hashTable, itr);
  }

  void increaseCount(packedKmerType kmer) {
    khint_t itr;
    int absent;
    itr = kc_c1_put(this->hashTable, kmer, &absent);
    if (absent) kh_val(this->hashTable, itr) = 0;
    ++kh_val(this->hashTable, itr);
  }

};

#endif //HASHTABLE_H
