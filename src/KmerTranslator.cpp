#include "KmerTranslator.h"

KmerTranslator::KmerTranslator(unsigned short span, unsigned short weight)
{
  this->span = span;
  this->weight = weight;

  //TODO: blockwise mask?
  //TODO: generalize!
  //11010111011011011001110011011011011101011
  this->_mask = (spacedKmerType) 1850513929963;
  this->_maskArray = new unsigned char[weight] {0,1,3,5,6,7,9,10,12,13,15,16,19,20,21,24,25,27,28,30,31,33,34,35,37,39,40};

    //110101110110110111011011011101011
    //this->_maskArray = new unsigned char[weight] {0,1,3,5,6,7,9,10,12,13,15,16,17,19,20,22,23,25,26,27,29,31,32};

  //this->_maskArray = new unsigned char[weight] {0,1,2,7,8,9,14,15,16,21,22,23,28,29,30,35,36,37,42,43,44,49,50,51,56,57,58};
    //this->_maskArray = new unsigned char[weight] {0,1,2,3,9,10,11,17,18,19,20,26,27,28,34,35,36,37,43,44,45,51,52,53,54};
}

unsigned short KmerTranslator::getSpan() const
{
  return this->span;
}

unsigned short KmerTranslator::getWeight() const
{
  return this->weight;
}

kmerType KmerTranslator::kmer2packedKmer(const spacedKmerType kmer) const
{
  kmerType packedKmer = 0;

  size_t j = 0;

  for (size_t i = 0; i < span; i++)
  {
     if (i == _maskArray[j])
     {
       packedKmer = (packedKmer<<2) | (kmerType) ((kmer & ((spacedKmerType)3 << (2*(span-1-i)))) >> (2*(span-1-i)));
       j++;
     }
   }
   
   return packedKmer;
}

kmerType KmerTranslator::kmer2minPackedKmer(const spacedKmerType kmer) const
{
  kmerType packedKmer = kmer2packedKmer(kmer);
  return(minIndex(packedKmer, weight));
}
