#include "uint128.h"

uint128_t largeInt2uint128(const LargeInt<1> value)
{
  return (uint128_t) value.getVal();
}

uint128_t largeInt2uint128(const LargeInt<2> value)
{
  return ((uint128_t)((value >> 64).getVal()) << 64) | (value.getVal());
}
