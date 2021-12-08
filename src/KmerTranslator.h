// Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>
#ifndef KMERTRANSLATOR_H
#define KMERTRANSLATOR_H

#include "kmer.h"

class KmerTranslator {

private:
  unsigned short span;   /* kmer size */
  unsigned short weight; /* number of informative positions */

  spacedKmerType _spaced_mask;
  unsigned char *_mask_array;
  unsigned char *_inverse_mask_array;

public:
  KmerTranslator(std::string spacedKmerPattern);

  ~KmerTranslator();

  packedKmerType kmer2packedKmer(const spacedKmerType kmer) const;

  packedKmerType kmer2minPackedKmer(const spacedKmerType kmer) const;

  packedKmerType kmer2minPackedKmer(const packedKmerType kmer) const;

  unsigned short getSpan() const;

  unsigned short getWeight() const;

  unsigned short getLongestBlock() const;

  void getBestSplit(unsigned int &logIndexSize, unsigned int &logOffsetSize) const;

  friend class CountProfile;
  
};

#endif // KMERTRANSLATOR_H
