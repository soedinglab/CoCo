//
// Created by aseidel1 on 7/21/20.
//
#include "LookuptableBase.h"
#include "khashl.h" // hash table

#ifndef COCO_HASHLOOKUPTABLE_H
#define COCO_HASHLOOKUPTABLE_H

KHASHL_MAP_INIT(, kc_c1_t, kc_c1, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)
class HashLookuptable: public LookupTableBase{

public:
    HashLookuptable(){
        hashTable = kc_c1_init();
    }

    kc_c1_t *hashTable;


    unsigned int getCount (const kmerType kmer) const {
        khint_t itr;
        int absent;
        itr = kc_c1_put(this->hashTable, kmer, &absent);
        if (absent) return 1;
        return kh_val(this->hashTable, itr);
    }

    void increaseCount(kmerType kmer){
        khint_t itr;
        int absent;
        itr = kc_c1_put(this->hashTable, kmer, &absent);
        if (absent) kh_val(this->hashTable, itr) = 0;
        ++kh_val(this->hashTable, itr);
    }



};
#endif //COCO_HASHLOOKUPTABLE_H
