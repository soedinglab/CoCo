// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#ifndef COUNTPROFILE_H
#define COUNTPROFILE_H

#include <assert.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "SequenceInfo.h"
#include "LookuptableBase.h"
#include "KmerTranslator.h"
#include "types.h"


struct __attribute__((__packed__)) CountProfileEntry {
  uint8_t valid;
  uint32_t count;          /* max count of spaced k-mers starting (matching) <readPos> */
};

class CountProfile {
private:

  CountProfileEntry *profile;
  const KmerTranslator *translator;
  const LookupTableBase *lookuptable;

  SequenceInfo *seqinfo;
  //const char* seqNuc;
  //const char* seqName;
  size_t maxprofileLength = 0,    /* maximal length of profile,
                                     corresponds to the allocated size of <profile> */
            profileLength = 0;    /* length of current profile,
                                     corresponds to the used size of <profile>*/



public:

  /* constructor */
  CountProfile(const KmerTranslator *translator, const LookupTableBase *lookuptable);

  /* destructor */
  ~CountProfile();

  /* create and store count profiles */
  void fill(SequenceInfo *seq, size_t length);

  /* getter */
  const char *getSeqName() { return (seqinfo->name.c_str()); }

  unsigned int getProfileLen() { return (profileLength); }

  SequenceInfo *getSeqInfo() { return (seqinfo); }

  /* basic profile operations */
  uint32_t* maximize();

  void addCountPerPosition(std::vector<uint32_t> &summedCountProfile);

  /* calculate abundance estimation value as 67% quantile over count profile */
  unsigned int calc67quantile();

  unsigned int calcMedian();

  unsigned int calcMedian(std::vector<uint32_t> &positionsOfInterest);

  unsigned int calcXquantile(double quantile, std::vector<uint32_t> &positionsOfInterest);

  /* advanced profile operations */

  bool checkForTransitionDrops(unsigned int minCount);
  bool checkForTransitionDropsNew(unsigned int minCount);
  bool checkForTransitionDropsNew2(unsigned int minCount);
  char checkForSpuriousTransitionDrops(uint32_t *maxProfile, unsigned int dropLevelCriterion, bool maskOnlyDropEdges=true);

  /* show tab-based table of positions and counts */
  void showProfile(FILE *fp = stdout) const;

};

#endif // COUNTPROFILE_H
