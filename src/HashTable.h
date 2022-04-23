// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>
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

  ~HashTable() {
    kc_c1_destroy(hashTable);
  }

  unsigned int getCount(const packedKmerType kmer) const {
    khint_t itr;
    int absent;
    itr = kc_c1_put(this->hashTable, kmer, &absent);
    if (absent) return 1;
    return kh_val(this->hashTable, itr);
  }

  bool increaseCount(packedKmerType kmer) {
    khint_t itr;
    int absent;
    itr = kc_c1_put(this->hashTable, kmer, &absent);
    if (absent) kh_val(this->hashTable, itr) = 0;
    ++kh_val(this->hashTable, itr);
    return true;
  }

  bool decreaseCount(packedKmerType kmer) {
    khint_t itr;
    int absent;
    itr = kc_c1_put(this->hashTable, kmer, &absent);
    if (absent) return false;
    --kh_val(this->hashTable, itr);
    return true;
  }

  void iterateOverAll(FILE* fp) const {
    EXIT(EXIT_FAILURE);
  }
  void printSize(){};
};

#endif //HASHTABLE_H
