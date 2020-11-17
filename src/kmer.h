// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#ifndef KMER
#define KMER

#include "uint128.h"

typedef uint128_t spacedKmerType;
typedef uint64_t packedKmerType;

packedKmerType minIndex(packedKmerType kmer, unsigned short kmerSize);

#endif // KMER


