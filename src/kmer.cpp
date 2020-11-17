// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include "kmer.h"
#include "emmintrin.h"
#include "immintrin.h"
#include "types.h"

#ifdef __SSSE3__


/* calculate reverse copmplement with SSE3 instructions
   Note: do not use typdef for packedKmerType here, because this function only works for uint64_t kmer
   if packedKmerType definition change use another function */
uint64_t revComplement_intr(uint64_t kmer, unsigned short k) {

  // broadcast 64bit to 128 bit
  __m128i x = _mm_cvtsi64_si128(kmer);

  // create lookup (set 16 bytes in 128 bit)
  // a lookup entry at the index of two nucleotids (4 bit) describes the reverse
  // complement of these two nucleotid in the higher 4 bits (lookup1) or in the
  // lower 4 bits (lookup2)
  __m128i lookup1 = _mm_set_epi8(0x50, 0x10, 0xD0, 0x90, 0x40, 0x00, 0xC0, 0x80, 0x70,
                                 0x30, 0xF0, 0xB0, 0x60, 0x20, 0xE0, 0xA0);
  __m128i lookup2 = _mm_set_epi8(0x05, 0x01, 0x0D, 0x09, 0x04, 0x00, 0x0C, 0x08, 0x07,
                                 0x03, 0x0F, 0x0B, 0x06, 0x02, 0x0E, 0x0A);

  // _mm_set1_epi8: create 128 bit with all bytes set to given value
  // here: 0x0F (00001111) and 0xF0 (11110000)
  // _mm_and_si128: bitwise AND
  __m128i kmer1 = _mm_and_si128(x, _mm_set1_epi8(0x0F)); // get lower 4 bits
  __m128i kmer2 = _mm_and_si128(x, _mm_set1_epi8(0xF0)); // get higher 4 bits

  // shift right by 2 nucleotids
  kmer2 >>= 4;

  // use _mm_shuffle_epi8 to look up reverse complement
  kmer1 = _mm_shuffle_epi8(lookup1, kmer1);
  kmer2 = _mm_shuffle_epi8(lookup2, kmer2);

  // _mm_or_si128: bitwise OR
  x = _mm_or_si128(kmer1, kmer2);

  // set upper 8 bytes to 0 and revert order of lower 8 bytes
  x = _mm_shuffle_epi8(x,
                       _mm_set_epi8(0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0, 1, 2, 3, 4, 5, 6, 7));

  // shift out the unused nucleotide positions (1 <= k <=32 )
  // broadcast 128 bit to 64 bit
  return (((uint64_t) _mm_cvtsi128_si64(x)) >> (uint64_t) (64 - 2 * k));
}

#endif

packedKmerType revComplement(packedKmerType kmer, unsigned short kmerSize) {

  packedKmerType revCompKmer = 0;

#ifdef __SSSE3__
  revCompKmer = revComplement_intr(kmer, kmerSize);
#else
  kmerType kmer_cp = kmer;
  kmerType mask = (packedKmerType) 3;

  for (int idx = 0; idx < kmerSize; idx++){
    kmerType nuc = (kmer_cp & mask);
    revCompKmer <<= 2;
    revCompKmer += int2rev[nuc];
    kmer_cp >>= 2;
  }
#endif

  return revCompKmer;
}


packedKmerType minIndex(packedKmerType kmer, unsigned short kmerSize) {
  packedKmerType revCompKmer = revComplement(kmer, kmerSize);

  return kmer < revCompKmer ? kmer : revCompKmer;
}
