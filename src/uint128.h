#ifndef UINT128
#define UINT128

#include <stdint.h>
#include <gatb/gatb_core.hpp>

#include <gatb/tools/math/LargeInt.hpp>

typedef __uint128_t uint128_t;

uint128_t largeInt2uint128(const LargeInt<1> value);

uint128_t largeInt2uint128(const LargeInt<2> value);

#endif // UINT128

