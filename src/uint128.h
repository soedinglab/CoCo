//Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#ifndef UINT128
#define UINT128

#include <stdint.h>
#include <gatb/gatb_core.hpp>

#include <gatb/tools/math/LargeInt.hpp>

typedef __uint128_t uint128_t;

uint128_t largeInt2uint128(const LargeInt<1> value);

uint128_t largeInt2uint128(const LargeInt<2> value);

std::ostream& operator<<( std::ostream& dest, __uint128_t value );

#endif // UINT128

