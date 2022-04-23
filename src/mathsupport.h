// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>

#ifndef MATHSUPPORT_H
#define MATHSUPPORT_H

#include <cstdlib>

inline size_t ipow(size_t base, size_t exponent) {
  size_t res = 1;
  if (base == 2) {
    return ((size_t) 1) << exponent;
  }

  for (size_t i = 0; i < exponent; i++)
    res = res * base;
  return res;
}

#endif // MATHSUPPORT_H

