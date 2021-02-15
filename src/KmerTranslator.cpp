// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include "KmerTranslator.h"

KmerTranslator::KmerTranslator() {
  span = 41;
  weight = 27;

  _span_mask = 2199023255551; //2^41-1
  //11010111011011011001110011011011011101011
  _mask_array = new unsigned char[weight]{0, 1, 3, 5, 6, 7, 9, 10, 12, 13, 15, 16, 19, 20, 21, 24, 25, 27, 28, 30, 31,
                                          33, 34, 35, 37, 39, 40};
  _inverse_mask_array = new unsigned char[span];
  unsigned short jdx = 0;
  for(unsigned short idx = 0; idx < span; idx++){
    if (_mask_array[jdx] ==idx) {
      _inverse_mask_array[idx] = jdx;
      jdx++;
    }
    else{
      _inverse_mask_array[idx] = UCHAR_MAX;
    }
  }
}

KmerTranslator::~KmerTranslator() {
  delete[] _mask_array;
}

unsigned short KmerTranslator::getSpan() const {
  return span;
}

unsigned short KmerTranslator::getWeight() const {
  return weight;
}

packedKmerType KmerTranslator::kmer2packedKmer(const spacedKmerType kmer) const {

  packedKmerType packedKmer = 0;

  size_t j = 0;

  for (size_t i = 0; i < span; i++) {
    if (i == _mask_array[j]) {
      packedKmer =
        (packedKmer << 2) | (packedKmerType) ((kmer & ((spacedKmerType) 3 << (2 * (span - 1 - i)))) >> (2 * (span - 1 - i)));
      j++;
    }
  }

  return packedKmer;
}

packedKmerType KmerTranslator::kmer2minPackedKmer(const spacedKmerType kmer) const {
  packedKmerType packedKmer = kmer2packedKmer(kmer);
  return (minIndex(packedKmer, weight));
}

packedKmerType KmerTranslator::kmer2minPackedKmer(const packedKmerType kmer) const {
  return (minIndex(kmer, weight));
}