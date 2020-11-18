// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#ifndef KMERTRANSLATOR_H
#define KMERTRANSLATOR_H

#include "kmer.h"

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

  friend class CountProfile;
  
};

#endif // KMERTRANSLATOR_H
