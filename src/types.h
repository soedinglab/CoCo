// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>

#ifndef TYPES
#define TYPES

#include <string>

typedef std::string SeqType;

// one character representation: A < C < T < G
const unsigned int alphabetSize = 4;
extern char int2res[alphabetSize];
extern char int2rev[alphabetSize];
extern int res2int['z' + 1];

void initialize();

#endif // TYPES

