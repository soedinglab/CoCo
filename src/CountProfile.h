#ifndef COUNTPROFILE_H
#define COUNTPROFILE_H

// Written by Annika Seidel

#include <assert.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "Lookuptable.h"
#include "KmerTranslator.h"
#include "types.h"



struct __attribute__((__packed__)) CountProfileEntry {
  uint8_t  nuc : 3;        /* nucleotid in 2 bit representation */
  uint16_t readPos : 10 ;  /* position in read */
  uint32_t count;          /* max count of spaced k-mers matching <readPos> */
};

class CountProfile
{
  private:

    CountProfileEntry*     profile;
    const KmerTranslator*  translator;
    const Lookuptable*     lookuptable;

    char *                 readName;             /* index of read */
    size_t                 maxprofileLength = 0, /* maximal length of profile,
                                                    corresponds to the allocated
                                                    size of <profile> */
                           profileLength = 0,    /* length of current profile,
                                                    corresponds to the used
                                                    size of <profile>*/
                           populationCoverage = 0;

public:

    /* constructor */
    CountProfile(const KmerTranslator *translator, const Lookuptable *lookuptable);

    CountProfile(const KmerTranslator *translator, const Lookuptable *lookuptable,
                 const SeqType seq, char *readName);

    /* destructor */
    ~CountProfile();

    /* create and store count profiles of <seq> */
    void fill(const SeqType seq, const char *readName);

    /* set estimate of population coverage c_pop to 67% quantile  */
    size_t calcPopulationCoverage();

    /* show tab-based table of readPos and abundance.
       abundance means here the maximal abundance of all spaced k-mers,
       which overlap readPos with a match */
    void showProfile(std::ostream &os=std::cout) const;

   //TODO: function: clear, correct
};

#endif // COUNTPROFILE_H
