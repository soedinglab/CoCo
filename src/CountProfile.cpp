// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include "CountProfile.h"
#include "mathsupport.h"
#include "Info.h"
#include <cstring>

#define WINDOW_SIZE 4


CountProfile::CountProfile(const KmerTranslator *translator,
                           const LookupTableBase *lookuptable) {
  this->translator = translator;
  this->lookuptable = lookuptable;
  this->profile = NULL;
  this->profile_length = 0;
  this->profile_length_alloc;
}

CountProfile::~CountProfile() {
  delete[] this->profile;
}

void CountProfile::fill(SequenceInfo *seqinfo, size_t length) {
  assert(this->lookuptable != NULL);
  assert(this->translator != NULL);

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  assert(length >= kmerSpan);
  this->seqinfo = seqinfo;
  this->profile_length = length - kmerSpan + 1;


  // check size of array and create new one in case of too small size
  if (profile_length > profile_length_alloc) {

    delete[] this->profile;
    this->profile = new CountProfileEntry[profile_length];
    this->profile_length_alloc = profile_length;
  }

  spacedKmerType kmer = 0, nStore = 0;
  for (size_t idx = 0; idx < length; idx++) {
    if ((char) res2int[(int) seqinfo->seq[idx]] != -1) {
      kmer = (kmer << 2) | (char) res2int[(int) seqinfo->seq[idx]];
      nStore = nStore << 1;
    } else // found non valid nucleotide
    {
      kmer = kmer << 2;
      nStore = nStore << 1 | 1;
    }

    if (idx >= kmerSpan - 1) {
      // check if k-mer contains 'N' (or any other invalid character)
      if  ((nStore & translator->_span_mask) != 0) {
        profile[idx - (kmerSpan - 1)].count = 0;
        profile[idx - (kmerSpan - 1)].valid = 0;
        continue;
      }

      packedKmerType packedKmer = translator->kmer2minPackedKmer(kmer);
      uint32_t count = lookuptable->getCount(packedKmer);
      profile[idx - (kmerSpan - 1)].count = count;
      profile[idx - (kmerSpan - 1)].valid = 1;
    }
  }


  //TODO: continue until maxprofillength to reset count and valid flag? not necessary
}


void CountProfile::showProfile(FILE *fp) const {
  for (size_t idx = 0; idx < profile_length; idx++) {
    fprintf(fp, "%u\t%u\n", idx, profile[idx].count);
  }
}

void CountProfile::showMaximzedProfile(FILE *fp) const {

  uint32_t *maxProfile = maximize();
  unsigned short kmerSpan = translator->getSpan();

  for (size_t idx = 0; idx < profile_length + kmerSpan - 1; idx++) {
    fprintf(fp, "%u\t%u\n", idx, maxProfile[idx]);
  }
}

uint32_t *CountProfile::maximize() const {

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  size_t maxProfileLen = profile_length + kmerSpan - 1;
  uint32_t *maxProfile = new uint32_t[maxProfileLen];
  for (size_t idx = 0; idx < maxProfileLen; idx++)
    maxProfile[idx] = 1;

  for (size_t idx = 0; idx < profile_length; idx++) {
    for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
      size_t pos = idx + translator->_mask_array[jdx];
      if (this->profile[idx].valid)
        maxProfile[pos] = std::max(this->profile[idx].count, maxProfile[pos]);
    }
  }

  return maxProfile;
}

void CountProfile::addCountPerPosition(std::vector<uint64_t> &summedCountProfile)
{
  if (summedCountProfile.size() < profile_length)
    summedCountProfile.resize(profile_length,0);

  for (size_t idx = 0; idx < profile_length; idx++)
    summedCountProfile[idx] += profile[idx].count;
}


unsigned int CountProfile::calcXquantile(double quantile, const std::vector<uint32_t> &positionsOfInterest) {

  unsigned int totalPos = 0 ;

  //copy count values
  uint32_t *counts;
  if (positionsOfInterest.empty()){
    totalPos = profile_length;
    counts = new uint32_t[totalPos];
    for (size_t idx = 0; idx < totalPos; idx++) {
      counts[idx] = profile[idx].count;
    }
  }
  else {
    for (size_t idx = 0; idx < positionsOfInterest.size(); idx++) {
      if (positionsOfInterest[idx] < profile_length)
        totalPos++;
    }
    counts = new uint32_t[totalPos];
    for (size_t idx = 0; idx < totalPos; idx++) {
      counts[idx] = profile[positionsOfInterest[idx]].count;
    }
  }

  // sort in ascending order
  sort(counts, counts + totalPos);

  auto x_quantile = counts[(uint32_t) (quantile * (double) totalPos)];

  delete[] counts;
  return x_quantile;
}

unsigned int CountProfile::calcMedian() {
  return(calcXquantile(0.5));
}

unsigned int CountProfile::calcMedian(const std::vector<uint32_t> &positionsOfInterest) {

  return(calcXquantile(0.5, positionsOfInterest));
}

bool CountProfile::correction(uint32_t *maxProfile, unsigned int covEst,  bool dryRun){

    unsigned short kmerSpan = this->translator->getSpan();
    unsigned short kmerWeight = translator->getWeight();
  unsigned int maxProfileLength = profile_length + kmerSpan - 1;
    unsigned int corrValues[maxProfileLength];
    memset(corrValues, 0, sizeof(*corrValues) * maxProfileLength);

    // correct values with probability for observing the same sequencing error multiple times
    std::vector<unsigned int> corrWindow;
    double corrFactor = 0.001;

    for(unsigned int idx = 0; idx < kmerSpan/2+1; idx++)
        corrWindow.push_back(maxProfile[idx]);
    corrValues[0] = corrFactor * *(std::max_element(corrWindow.begin(), corrWindow.end())) + 1;

    for(unsigned int idx = 1; idx < maxProfileLength; idx++) {
        if (idx + kmerSpan/2 < maxProfileLength)
            corrWindow.push_back(maxProfile[idx + kmerSpan/2]);
        if (idx > kmerSpan/2)
            corrWindow.erase(corrWindow.begin());
        unsigned int max = *(std::max_element(corrWindow.begin(), corrWindow.end()));
        corrValues[idx] = (unsigned int)(corrFactor * max + 1);

    }

    // 1. find sequencing errors

    unsigned int candidates[maxProfileLength];
    bool foundErrorPos = false;
    memset(candidates, 0, sizeof(*candidates) * maxProfileLength);

    for (size_t idx = 0; idx < maxProfileLength; idx++){
        if (maxProfile[idx] <= 1 + corrValues[idx] && maxProfile[idx] <= 0.1 * covEst) {
            candidates[idx] = 1;
            if (dryRun && maxProfile[idx] <= 0.1 * covEst) {//TODO: change 10% later for correction
                 foundErrorPos = true;
                 indexVec.push_back(idx);
                 Info(Info::DEBUG) << seqinfo->name.c_str() << "\t" << idx << "\n";
            }
        }
    }

    //if (dryRun)
    return foundErrorPos;
    

    //  2. correct sequencing errors
    // TODO: 1. majority voting? 2. local covEst 3. global covEst (= median)

    //return false;
}

bool CountProfile::checkForSpuriousTransitionDrops(uint32_t *maxProfile, unsigned int dropLevelCriterion, bool maskOnlyDropEdges) {

  uint8_t candidates[profile_length];
  memset(candidates, 0, sizeof(*candidates) * profile_length);

  double correctionFactor= 0.001;
  unsigned int dropstart = 0;
  unsigned int dropend = profile_length;

  // 1) find drops -> candidates

  for (size_t idx = 1; idx < profile_length; idx++) {

    // dropstart
    if (this->profile[idx].count / this->profile[idx - 1].count <= 0.5 &&
        this->profile[idx].count <= dropLevelCriterion &&
        this->profile[idx - 1].count > dropLevelCriterion) {
      dropstart = idx;
    }
      // dropend
    else if (this->profile[idx - 1].count / this->profile[idx].count <= 0.5 &&
             this->profile[idx - 1].count <= dropLevelCriterion &&
             this->profile[idx].count > dropLevelCriterion) {

      if (dropstart == profile_length)
        continue;
      /*if (dropstart == 0)
        continue;*/

      dropend = idx;

      for (size_t d = dropstart; d < dropend; d++) {
        if (profile[d].count <= dropLevelCriterion)
          candidates[d] = 1;
      }

      dropstart = profile_length;
      dropend = profile_length;
    }
  }

  if(dropstart>0 && dropstart < profile_length) {

    for (size_t d = dropstart; d < profile_length; d++) {
      if (profile[d].count <= dropLevelCriterion)
        candidates[d] = 1;
    }
  }

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  size_t maxProfileLen = profile_length + kmerSpan - 1;
  dropstart = 0;
  dropend = maxProfileLen;

  /* 2) check if drop remains on maximized profile -> normal drop (spurious position)
        or if drop disappear on maximized profile -> transition drop (spurious order) */

  for (size_t idx = 1; idx < maxProfileLen; idx++) {

    // dropstart
    if (maxProfile[idx] / maxProfile[idx-1] < 0.5 &&
        maxProfile[idx] <= (dropLevelCriterion + (correctionFactor*maxProfile[idx-1]+1)) &&
        maxProfile[idx-1] > (dropLevelCriterion + (correctionFactor*maxProfile[idx-1]+1))) {
      dropstart = idx;
    }
      // dropend
    else if (maxProfile[idx-1] / maxProfile[idx] < 0.5 &&
             maxProfile[idx-1] <= (dropLevelCriterion + (correctionFactor*maxProfile[idx]+1)) &&
             maxProfile[idx] > (dropLevelCriterion + (correctionFactor*maxProfile[idx]+1))) {

      if (dropstart == maxProfileLen)
        continue;

      dropend = idx;
      if (maskOnlyDropEdges) {

        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = dropstart - translator->_mask_array[jdx];
          if (pos >= 0 && pos < profile_length)
            candidates[pos] = 0;
        }

        if (dropend > dropstart + 1) {
          for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            int pos = dropend - 1 - translator->_mask_array[jdx];
            if (pos >= 0 && pos < profile_length)
              candidates[pos] = 0;
          }
        }
      }
      else{
        for (size_t d = dropstart; d < dropend; d++) {
          for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            int pos = d - translator->_mask_array[jdx];
            if (pos >= 0 && pos < profile_length)
              candidates[pos] = 0;
          }
        }
      }
      dropstart=maxProfileLen;
      dropend=maxProfileLen;
    }
  }

  if(dropstart>0 && dropstart < maxProfileLen) {

    if (maskOnlyDropEdges) {
      for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
        int pos = dropstart - translator->_mask_array[jdx];
        if (pos >= 0 && pos < profile_length)
          candidates[pos] = 0;
      }

      if (dropend > dropstart + 1) {
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = dropend - 1 - translator->_mask_array[jdx];
          if (pos >= 0 && pos < profile_length)
            candidates[pos] = 0;
        }
      }
    }
    else {
      for (size_t d = dropstart; d < maxProfileLen; d++) {
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = d - translator->_mask_array[jdx];
          if (pos >= 0 && pos < profile_length)
            candidates[pos] = 0;
        }
      }
    }
  }

  for (size_t idx = 0; idx < this->profile_length; idx++) {
    if (candidates[idx] != 0)
      return true;
  }

  return false;

  //TODO:N's
}


bool CountProfile::checkForSpuriousTransitionDropsWithWindow(uint32_t *maxProfile, unsigned int covEst,
                                                             double localPercDrop, double globalPercDrop,
                                                             bool maskOnlyDropEdges)
{

  unsigned short kmerSpan = this->translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  unsigned int corrValues[profile_length];
  memset(corrValues, 0, sizeof(*corrValues) * profile_length);

  // correct values with probability for observing the same sequencing error multiple times
  std::vector<unsigned int> corrWindow;
  double corrFactor = 0.001;

  for(unsigned int idx = 0; idx < kmerSpan/2+1; idx++)
    corrWindow.push_back(this->profile[idx].count);
  corrValues[0] = corrFactor * *(std::max_element(corrWindow.begin(), corrWindow.end())) + 1;

  for(unsigned int idx = 1; idx < profile_length; idx++) {
    if (idx + kmerSpan/2 < profile_length)
      corrWindow.push_back(this->profile[idx + kmerSpan/2].count);
    if (idx > kmerSpan/2)
      corrWindow.erase(corrWindow.begin());
    unsigned int max = *(std::max_element(corrWindow.begin(), corrWindow.end()));
    corrValues[idx] = (unsigned int)(corrFactor * max+1);

  }

  unsigned int candidates[profile_length];
  memset(candidates, 0, sizeof(*candidates) * profile_length);

  unsigned int maxProfileLength = profile_length + kmerSpan - 1;

  vector<uint32_t> window_back;
  vector<uint32_t> window_front;

  // fill windows for first entry
  window_back.push_back(this->profile[0].count);
  for (size_t idx = 1; idx < WINDOW_SIZE+1; idx++){
    window_front.push_back(this->profile[idx].count);
  }

  uint32_t dropstart_level = UINT_MAX, dropend_level = UINT_MAX, compare_level;
  unsigned int dropstart = profile_length, dropend = profile_length;

  /* 1) find candidates */

  if(this->profile[0].count < covEst * globalPercDrop)
    dropstart = 0;

  for (size_t idx = 1; idx < profile_length; idx++){

    uint32_t window_back_level = *(std::min_element(std::begin(window_back), std::end(window_back)));
    uint32_t window_front_level = *(std::min_element(std::begin(window_front), std::end(window_front)));

    if ((double) this->profile[idx].count < localPercDrop * ((double) window_back_level-(double)corrValues[idx-1]) &&\
        window_back_level - this->profile[idx].count >= SIGNIFICANT_LEVEL_DIFF && \
        this->profile[idx].count < (double)covEst * globalPercDrop){

      //drop start
      dropstart = idx;
      dropend = profile_length;
      dropstart_level = window_back_level;
    }

    else if (window_front.size() == WINDOW_SIZE && \
            ((localPercDrop * ((double) window_front_level - (double) corrValues[idx] ) >  (double) window_back_level
              && window_front_level - window_back_level >= SIGNIFICANT_LEVEL_DIFF) || window_front_level > (double)covEst * globalPercDrop)){

      //drop end (only considered if we already observed a dropstart or if dropstart=0)
      if (dropstart != this->profile_length) {
        dropend = idx;

        // mark candidates
        dropend_level = window_front_level;
        compare_level = std::min(dropstart_level, dropend_level);
        // set all positions within drop boundaries as candidates if they are smaller than percDrop * compare_level
        for (size_t d = dropstart; d < dropend; d++) {
          if (this->profile[d].count < localPercDrop * compare_level) {
            candidates[d] = compare_level;
          }
        }

        /* 2) check if drop remains on maximized profile -> normal drop (spurious position)
              or if drop disappear on maximized profile -> transition drop (spurious order) */

        // check candidates
        for (size_t d = dropstart; d < dropend; d++) {
          if (maxProfile[d] < localPercDrop * compare_level && (!maskOnlyDropEdges || \
             ((d > 0 && maxProfile[d] < localPercDrop * maxProfile[d - 1]) || maxProfile[d] < localPercDrop * maxProfile[d + 1]))) {
            for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
              int pos = d - translator->_mask_array[jdx];
              if (pos >= dropstart && pos < profile_length)
                candidates[pos] = 0;
            }
          }
        }
      }

      // reset values
      dropstart = profile_length;
      dropstart_level = UINT_MAX;
      dropend_level = UINT_MAX;
      dropend = profile_length;
    }

    // update windows
    if (window_back.size()>= WINDOW_SIZE)
      window_back.erase(window_back.begin());
    window_back.push_back(this->profile[idx].count);

    window_front.erase(window_front.begin());
    if (idx + WINDOW_SIZE < profile_length)
      window_front.push_back(this->profile[idx+WINDOW_SIZE].count);
  }

  // edge case: check if there is a drop region until the end of the read (dropstart without dropend)
  if(dropstart != 0 && dropstart != profile_length){
    for (size_t d = dropstart; d < dropend; d++) {
      if (this->profile[d].count < localPercDrop * dropstart_level) {
        candidates[d] = dropstart_level;
      }
    }
    for (size_t d = dropstart; d < maxProfileLength; d++) {
      if (maxProfile[d] < localPercDrop * dropstart_level && (!maskOnlyDropEdges ||\
         ((d > 0 && maxProfile[d] < localPercDrop * maxProfile[d - 1]) || \
          (d + 1 < maxProfileLength &&  maxProfile[d] < localPercDrop * maxProfile[d + 1])))) {
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = d - translator->_mask_array[jdx];
          if (pos >= dropstart)
            if(pos < profile_length)
              candidates[pos] = 0;
        }
      }
    }
  }

  // check for unexplained drops -> transition drops
  for(unsigned int idx=0; idx < profile_length; idx++){
    if (candidates[idx] != 0)
      return true;
  }

  return false;
}

std::vector<int> CountProfile::getIdx() {
  return indexVec;
}

void CountProfile::resetIdx() {
  indexVec = {};
}

