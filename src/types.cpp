#include "types.h"

char int2res[alphabetSize];
char int2rev[alphabetSize];
int res2int['Z'+1];

void initialize()
{
  int2res[0] = 'A'; // T
  int2res[1] = 'C'; // G
  int2res[2] = 'T'; // A
  int2res[3] = 'G'; // C
  int2rev[0] = 2;
  int2rev[1] = 3;
  int2rev[2] = 0;
  int2rev[3] = 1;

  for (unsigned int i = 0; i <= 'Z'; ++i)
  {
    res2int[i]=-1;
  }
  for (unsigned int i = 0; i < alphabetSize; ++i)
  {
    res2int[(int)int2res[i]] = i;
  }
}
