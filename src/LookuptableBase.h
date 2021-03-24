// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include "types.h"
#include "kmer.h"

#ifndef LOOKUPTABLEBASE_H
#define LOOKUPTABLEBASE_H

class LookupTableBase {
public:
  virtual ~LookupTableBase() {};
  virtual unsigned int getCount(const packedKmerType kmer) const = 0;
  virtual void iterateOverAll(FILE* fp) const = 0;
  virtual bool decreaseCount(packedKmerType kmer) = 0;

  virtual bool increaseCount(packedKmerType kmer) = 0;
  virtual  void printSize() =0;
};

#endif //LOOKUPTABLEBASE_H
