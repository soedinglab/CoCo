//
// Created by aseidel1 on 7/21/20.
//
#include "types.h"
#include "kmer.h"
#ifndef COCO_LOOKUPTABLEBASE_H
#define COCO_LOOKUPTABLEBASE_H

class LookupTableBase{
public:
    virtual unsigned int getCount (const kmerType kmer) const = 0;
};

#endif //COCO_LOOKUPTABLEBASE_H
