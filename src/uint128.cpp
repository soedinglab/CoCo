//Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>

#include "uint128.h"

/*uint128_t largeInt2uint128(const LargeInt<1> value)
{
  return (uint128_t) value.getVal();
}*/

uint128_t largeInt2uint128(const LargeInt<2> value) {
  return ((uint128_t) ((value >> 64).getVal()) << 64) | (value.getVal());
}

//https://stackoverflow.com/questions/25114597/how-to-print-int128-in-g
std::ostream &operator<<(std::ostream &dest, __uint128_t value) {
  std::ostream::sentry s(dest);
  if (s) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[128];
    char *d = std::end(buffer);
    do {
      --d;
      *d = "0123456789"[tmp % 10];
      tmp /= 10;
    } while (tmp != 0);
    if (value < 0) {
      --d;
      *d = '-';
    }
    int len = std::end(buffer) - d;
    if (dest.rdbuf()->sputn(d, len) != len) {
      dest.setstate(std::ios_base::badbit);
    }
  }
  return dest;
}

