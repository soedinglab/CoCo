#ifndef COUNTPROFILE_H
#define COUNTPROFILE_H

// Written by Annika Seidel

#include <assert.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "LookuptableBase.h"
#include "KmerTranslator.h"
#include "types.h"



struct __attribute__((__packed__)) CountProfileEntry {
  uint8_t  nuc : 3;        /* nucleotid in 2 bit representation */
  uint8_t valid;
  uint32_t count;          /* max count of spaced k-mers matching <readPos> */
};

class CountProfile
{
  private:

    CountProfileEntry*     profile;
    const KmerTranslator*  translator;
    const LookupTableBase*     lookuptable;

    char *                 seqName;              /* index of sequence */
    size_t                 maxprofileLength = 0, /* maximal length of profile,
                                                    corresponds to the allocated
                                                    size of <profile> */
                           profileLength = 0,    /* length of current profile,
                                                    corresponds to the used
                                                    size of <profile>*/
                           abundanceEstimation = 0;

public:

    /* constructor */
    CountProfile(const KmerTranslator *translator, const LookupTableBase *lookuptable);

    CountProfile(const KmerTranslator *translator, const LookupTableBase *lookuptable,
                 const SeqType seq, char *seqdName);

    /* destructor */
    ~CountProfile();

    /* getter */
    char *getSeqName(){return (this->seqName);}
    float getCorrFactor(){return (0);}
    /* create and store count profiles of <seq> */
    void fill(const SeqType seq, const char *seqName);

    /* set abundanceEstimation value c_pop to 67% quantile over count profile */
    size_t calc67quantile();

    std::vector<unsigned int> getDropPointsInMaximzedProfile();

    bool  checkForRiseAndDropPoints (std::vector<unsigned int> dropPositions, unsigned int windowsize);


    /* show tab-based table of seqPos and count.
       count means here the maximal count of all spaced k-mers,
       which overlap seqPos with a match */
    void showProfile(FILE *fp=stdout) const;

   //TODO: function: clear, correct
    size_t calcMedian();
};

#endif // COUNTPROFILE_H
