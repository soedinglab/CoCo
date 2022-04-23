// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>
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
  packedKmerType kmer;
};

#define ERROR_FREE 0
#define NONE_CORRECTED 1
#define SOME_CORRECTED 2
#define ALL_CORRECTED 3
#define TOO_MANY_ERRORS 4

class CountProfile {
private:

  CountProfileEntry *profile;
  const KmerTranslator *translator;
  LookupTableBase *lookuptable;

  SequenceInfo *seqinfo;

  uint32_t profile_length_alloc = 0,  /* corresponds to the allocated size of <profile> */
           profile_length = 0;        /* length of current profile, corresponds to the used size of <profile>*/


public:

  int ambigCorr;

  /* constructor */
  CountProfile(const KmerTranslator *translator, LookupTableBase *lookuptable);

  /* destructor */
  ~CountProfile();

  /* create and store count profiles */
  void fill(SequenceInfo *seq);

  void update(bool updateLookupTable=false);

  /* setter */
  //set seqinfo, but do not update countprofile, mark this by profile_length=0
  void setSeqInfo(SequenceInfo *seqinfo);

  /* getter */
  const char *getSeqName() { return (seqinfo->name.c_str()); }

  unsigned int getProfileLen() { return (profile_length); }

  SequenceInfo *getSeqInfo() { return (seqinfo); }

  CountProfileEntry * getProfilecopy();

  CountProfileEntry * getProfile();

  /*** show functions ***/

  void showProfile(FILE *fp = stdout) const;

  void showMaximzedProfile(FILE *fp = stdout) const;

  /*** basic profile operations ***/
  uint32_t* maximize() const;

  void addCountPerPosition(std::vector<uint64_t> &summedCountProfile);

  unsigned int calcXquantile(double quantile, const std::vector<uint32_t> &positionsOfInterest = std::vector<uint32_t>());

  unsigned int calcMedian();

  unsigned int calcMedian(const std::vector<uint32_t> &positionsOfInterest);

  /*** advanced profile operations ***/

  /* correction operations */

  int doSubstitutionCorrection(uint32_t *maxProfile, double threshold, unsigned int pseudocount,
                               unsigned int lowerbound, bool needMultipleKmers, bool updateLookup, unsigned int *correctedSubstitutions);

  bool doIndelCorrection(uint32_t *maxProfile, double threshold, unsigned int pseudocount, unsigned int lowerbound, bool trySubstitution,
                         bool updateLookup, unsigned int *correctedSubstitutions, unsigned int *correctedInsertions, unsigned int *correctedDeletions);

  bool doTrimming(uint32_t *maxProfile, double threshold, unsigned int pseudocount, unsigned int lowerbound,
                  unsigned int maxTrimLen, bool updateLookup, unsigned int *trimmedCounter);

  int firstLastUniqueKmerCorrectionStrategy(unsigned int substitutionStart, unsigned int firstUniqueKmerStart,
                                            unsigned int lastUniqueKmerStart, const uint32_t *neighborhoodTolerance);

  int edgeSubstitutionCorrection(unsigned int errorPos, uint32_t *neighborhoodTolerance);

  bool tryInsertionCorrection(unsigned int insertionStart, unsigned int insertionLen, const uint32_t *neighborhoodTolerance);

  int tryDeletionCorrection(unsigned int deletionPos, const uint32_t *neighborhoodTolerance);

  /* filter operations */

  bool checkForSpuriousTransitionDropsWithWindowNew(double threshold);

  bool checkForSpuriousTransitionDropsWithWindow(uint32_t *maxProfile, unsigned int covEst, double localPercDrop, \
                                                 double globalPercDrop, bool maskOnlyDropEdges=true);
  //outdated
  bool checkForSpuriousTransitionDrops(uint32_t *maxProfile, unsigned int dropLevelCriterion, bool maskOnlyDropEdges=true);

};

#endif // COUNTPROFILE_H
