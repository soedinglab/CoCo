#ifndef KMERTRANSLATOR_H
#define KMERTRANSLATOR_H

#include "kmer.h"
class KmerTranslator
{

private:
  unsigned short span;   /* kmer size */
  unsigned short weight; /* number of informative positions */

  spacedKmerType _mask;
  unsigned char* _maskArray;

public:
  KmerTranslator(unsigned short span, unsigned short weight);

  kmerType kmer2packedKmer(const spacedKmerType kmer) const;
  kmerType kmer2minPackedKmer(const spacedKmerType kmer) const;
  unsigned short getSpan() const;
  unsigned short getWeight() const;
};

#endif // KMERTRANSLATOR_H
