// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>
#ifndef KMER
#define KMER

#include "uint128.h"

typedef uint128_t spacedKmerType;
typedef uint64_t packedKmerType;

packedKmerType minIndex(packedKmerType kmer, unsigned short kmerSize);

char* packedKmer2String(packedKmerType kmer, unsigned short kmerSize);


#endif // KMER


