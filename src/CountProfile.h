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

#define SIGNIFICANT_LEVEL_DIFF 10


struct __attribute__((__packed__)) CountProfileEntry {
  uint8_t valid;
  uint32_t count;
};

class CountProfile {
private:

  CountProfileEntry *profile;
  const KmerTranslator *translator;
  const LookupTableBase *lookuptable;

  SequenceInfo *seqinfo;

  uint32_t profile_length_alloc = 0,  /* corresponds to the allocated size of <profile> */
           profile_length = 0;        /* length of current profile, corresponds to the used size of <profile>*/

public:

  /* constructor */
  CountProfile(const KmerTranslator *translator, const LookupTableBase *lookuptable);

  /* destructor */
  ~CountProfile();

  /* create and store count profiles */
  void fill(SequenceInfo *seq, size_t length);

  /* getter */
  const char *getSeqName() { return (seqinfo->name.c_str()); }

  unsigned int getProfileLen() { return (profile_length); }

  SequenceInfo *getSeqInfo() { return (seqinfo); }

  /* show tab-based table of positions and counts */
  void showProfile(FILE *fp = stdout) const;
  void showMaximzedProfile(FILE *fp = stdout) const;

  /*** basic profile operations ***/
  uint32_t* maximize() const;

  void addCountPerPosition(std::vector<uint64_t> &summedCountProfile);

  unsigned int calcXquantile(double quantile, const std::vector<uint32_t> &positionsOfInterest = std::vector<uint32_t>());

  unsigned int calcMedian();

  unsigned int calcMedian(const std::vector<uint32_t> &positionsOfInterest);

  /*** advanced profile operations ***/

  bool correction(uint32_t *maxProfile, unsigned int covEst,  bool dryRun);

  bool checkForSpuriousTransitionDropsWithWindow(uint32_t *maxProfile, unsigned int covEst, double localPercDrop, \
                                                 double globalPercDrop, bool maskOnlyDropEdges=true);
  //outdated
  bool checkForSpuriousTransitionDrops(uint32_t *maxProfile, unsigned int dropLevelCriterion, bool maskOnlyDropEdges=true);


};

#endif // COUNTPROFILE_H
