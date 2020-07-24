// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
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

typedef struct{
  string  name, comment, seq, qual;
  char sep;
}
SequenceInfo;


struct __attribute__((__packed__)) CountProfileEntry {
  uint8_t valid;
  uint32_t count;          /* max count of spaced k-mers matching <readPos> */
};

class CountProfile
{
  private:

    CountProfileEntry* profile;
    const KmerTranslator* translator;
    const LookupTableBase* lookuptable;

    SequenceInfo* seqinfo;
    //const char* seqNuc;
    //const char* seqName;
    size_t  maxprofileLength = 0, /* maximal length of profile,
                                     corresponds to the allocated size of <profile> */
            profileLength = 0;    /* length of current profile,
                                     corresponds to the used size of <profile>*/

public:

    /* constructor */
    CountProfile(const KmerTranslator *translator, const LookupTableBase *lookuptable);

    /* destructor */
    ~CountProfile();

    /* getter */
    const char *getSeqName(){return (seqinfo->name.c_str());}
    SequenceInfo *getSeqInfo(){return (seqinfo);}

    /* create and store count profiles */
    void fill(SequenceInfo *seq, size_t length);

    /* calculate abundanceEstimation value as 67% quantile over count profile */
    unsigned int calc67quantile();

    std::vector<unsigned int> getDropPointsInMaximzedProfile();

    bool  checkForRiseAndDropPoints (std::vector<unsigned int> dropPositions);

    /* show tab-based table of seqPos and count */
    void showProfile(FILE *fp = stdout) const;

};

#endif // COUNTPROFILE_H
