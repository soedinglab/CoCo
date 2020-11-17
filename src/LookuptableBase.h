// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include "types.h"
#include "kmer.h"

#ifndef LOOKUPTABLEBASE_H
#define LOOKUPTABLEBASE_H

class LookupTableBase {
public:
  virtual unsigned int getCount(const packedKmerType kmer) const = 0;
};

#endif //LOOKUPTABLEBASE_H
