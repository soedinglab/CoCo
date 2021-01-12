// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#ifndef KMERTRANSLATOR_H
#define KMERTRANSLATOR_H

#include "kmer.h"

//Debug
#include <iostream>

class KmerTranslator {

private:
  unsigned short span;   /* kmer size */
  unsigned short weight; /* number of informative positions */

  spacedKmerType _span_mask;
  unsigned char *_mask_array;

public:
  KmerTranslator();

  ~KmerTranslator();

  packedKmerType kmer2packedKmer(const spacedKmerType kmer) const;

  packedKmerType kmer2minPackedKmer(const spacedKmerType kmer) const;

  unsigned short getSpan() const;

  unsigned short getWeight() const;

  void setSW(unsigned short newSpan, unsigned short newWeight);

  void setMask(unsigned char *newMask);

  friend class CountProfile;

  //debug
  void printmask(std::string prefix);
  
};

#endif // KMERTRANSLATOR_H
