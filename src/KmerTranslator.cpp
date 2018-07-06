#include "KmerTranslator.h"

KmerTranslator::KmerTranslator(unsigned short span, unsigned short weight)
{
  this->span = span;
  this->weight = weight;

  //TODO: generalize!
  //11010111011011011001110011011011011101011
  this->_mask = (spacedKmerType) 1850513929963;
  this->_maskArray = new unsigned char[weight] {0,1,3,5,6,7,9,10,12,13,15,16,19,20,21,24,25,27,28,30,31,33,34,35,37,39,40};
}

kmerType KmerTranslator::kmer2packedKmer(const spacedKmerType kmer) const
{
  kmerType packedKmer = 0;

  size_t j = 0;

  for (size_t i = 0; i < span; i++)
  {
     if (i == _maskArray[j])
     {
       packedKmer = (packedKmer<<2) | (kmerType) ((kmer & ((__uint128_t)3 << (2*(span-1-i)))) >> (2*(span-1-i)));
       j++;
     }
   }

   return packedKmer;
}

kmerType KmerTranslator::kmer2minPackedKmer(const spacedKmerType kmer) const
{
  kmerType packedKmer = kmer2packedKmer(kmer);
  return(minIndex(packedKmer, span));
}