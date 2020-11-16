// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include "CountProfile.h"
#include "mathsupport.h"
#include <cstring>

#define WINDOW_SIZE 4
#define SIGNIFICANT_LEVEL_DIFF 10

#define MIN_UPPER_LEVEL_POSITIONS 3
#define MIN_UPPER_LEVEL_POSITIONS_NEW 1

CountProfile::CountProfile(const KmerTranslator *translator,
                           const LookupTableBase *lookuptable) {
  this->translator = translator;
  this->lookuptable = lookuptable;
  this->profile = NULL;
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
  this->profileLength = length - kmerSpan + 1;


  // check size of array and create new one in case of too small size
  if (this->profileLength > this->maxprofileLength) {
    delete[] this->profile;
    this->profile = new CountProfileEntry[this->profileLength];
    this->maxprofileLength = this->profileLength;
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
      // check if spaced k-mer contains 'N'
      if ((nStore & (spacedKmerType) 2199023255551) != 0) { //if ((nStore & translator->_mask) != 0)
        profile[idx - (kmerSpan - 1)].count = 0;
        profile[idx - (kmerSpan - 1)].valid = 0;
        continue;
      }

      kmerType packedKmer = translator->kmer2minPackedKmer(kmer);
      uint32_t count = lookuptable->getCount(packedKmer);
      profile[idx - (kmerSpan - 1)].count = count;
      profile[idx - (kmerSpan - 1)].valid = 1;
    }
  }

  /*std::cout << "count profile" << std::endl;
  for (size_t idx = 0; idx < this->profileLength; idx++) {
    std::cout << idx << ":" << this->profile[idx].count << std::endl;
  }*/
  //TODO: continue until maxprofillength to reset count and valid flag? not necessary
}

uint32_t *CountProfile::maximize() {

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  size_t maxProfileLen = this->profileLength + kmerSpan - 1;
  uint32_t *maxProfile = new uint32_t[maxProfileLen];
  for (size_t idx = 0; idx < maxProfileLen; idx++)
    maxProfile[idx] = 1;

  for (size_t idx = 0; idx < this->profileLength; idx++) {
    for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
      size_t pos = idx + translator->_maskArray[jdx];
      if (this->profile[idx].valid)
        maxProfile[pos] = std::max(this->profile[idx].count, maxProfile[pos]);
    }
  }

  /*std::cout << "maximized count profile" << std::endl;
  for (size_t idx = 0; idx < maxProfileLen; idx++) {
   std::cout << idx << ":" << maxProfile[idx] << std::endl;
  }*/

  return maxProfile;
}

void CountProfile::addCountPerPosition(std::vector<uint32_t> &summedCountProfile)
{
  if (summedCountProfile.size() < this->profileLength)
    summedCountProfile.resize(this->profileLength,0);

  for (size_t idx = 0; idx < this->profileLength; idx++)
   summedCountProfile[idx] += profile[idx].count;
}



void CountProfile::showProfile(FILE *fp) const {
  for (size_t idx = 0; idx < profileLength; idx++) {
    fprintf(fp, "%u\t%u\n", idx, profile[idx].count);
  }

  /*unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  size_t maxProfileLen = this->profileLength + kmerSpan - 1;
  uint32_t maxProfile[maxProfileLen];
  for (size_t idx = 0; idx < maxProfileLen; idx++)
    maxProfile[idx] = 1;

  for (size_t idx = 0; idx < this->profileLength; idx++) {
    for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
      size_t pos = idx + translator->_maskArray[jdx];
      if (this->profile[idx].valid)
        maxProfile[pos] = std::max(this->profile[idx].count, maxProfile[pos]);
    }
  }
  for (size_t idx = 0; idx < maxProfileLen; idx++) {
    fprintf(fp, "%u\t%u\n", idx, maxProfile[idx]);
  }*/
}


unsigned int CountProfile::calcXquantile(double quantile, std::vector<uint32_t> &positionsOfInterest) {
  unsigned int totalPos = 0 ;
  for (size_t idx = 0; idx < positionsOfInterest.size(); idx++) {
    if (positionsOfInterest[idx] < profileLength)
      totalPos++;
  }
  //copy count values
  uint32_t *counts = new uint32_t[totalPos];

  for (size_t idx = 0; idx < totalPos; idx++) {
    counts[idx] = profile[positionsOfInterest[idx]].count;
  }

  // sort in ascending order
  sort(counts, counts + totalPos);

  auto x_quantile = counts[(uint32_t) (quantile * (double) totalPos)];

  delete[] counts;
  return x_quantile;
}

unsigned int CountProfile::calcMedian() {
  //copy count values
  uint32_t *counts = new uint32_t[profileLength];

  for (size_t idx = 0; idx < profileLength; idx++) {
    counts[idx] = profile[idx].count;
  }

  // sort in ascending order
  sort(counts, counts + profileLength);

  auto median = counts[(uint32_t) (0.5 * (double) this->profileLength)];

  delete[] counts;
  return median;
}

unsigned int CountProfile::calcMedian(std::vector<uint32_t> &positionsOfInterest) {

  return(calcXquantile(0.5, positionsOfInterest));
}


char CountProfile::checkForSpuriousTransitionDrops(uint32_t *maxProfile, unsigned int dropLevelCriterion, bool maskOnlyDropEdges) {

  uint8_t candidates[this->profileLength];
  memset(candidates, 0, sizeof(*candidates) * this->profileLength);

  double correctionFactor= 0.01;
  unsigned int dropstart = 0;
  unsigned int dropend = this->profileLength;
  for (size_t idx = 1; idx < this->profileLength; idx++) {

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

      if (dropstart == this->profileLength)
        continue;
      /*if (dropstart == 0)
        continue;*/

      dropend = idx;

      for (size_t d = dropstart; d < dropend; d++) {
        if (profile[d].count <= dropLevelCriterion)
          candidates[d] = 1;
      }

      dropstart=this->profileLength;
      dropend = this->profileLength;
    }
  }

  if(dropstart>0 && dropstart < this->profileLength) {

    for (size_t d = dropstart; d < this->profileLength; d++) {
      if (profile[d].count <= dropLevelCriterion)
        candidates[d] = 1;
    }
  }

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  size_t maxProfileLen = this->profileLength + kmerSpan - 1;

  dropstart = 0;
  dropend = maxProfileLen;

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
          int pos = dropstart - translator->_maskArray[jdx];
          if (pos >= 0 && pos < profileLength)
            candidates[pos] = 0;
        }

        if (dropend > dropstart + 1) {
          for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            int pos = dropend - 1 - translator->_maskArray[jdx];
            if (pos >= 0 && pos < profileLength)
              candidates[pos] = 0;
          }
        }
      }
      else{
        for (size_t d = dropstart; d < dropend; d++) {
          for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            int pos = d - translator->_maskArray[jdx];
            if (pos >= 0 && pos < profileLength)
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
        int pos = dropstart - translator->_maskArray[jdx];
        if (pos >= 0 && pos < profileLength)
          candidates[pos] = 0;
      }

      if (dropend > dropstart + 1) {
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = dropend - 1 - translator->_maskArray[jdx];
          if (pos >= 0 && pos < profileLength)
            candidates[pos] = 0;
        }
      }
    }
    else {
      for (size_t d = dropstart; d < maxProfileLen; d++) {
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = d - translator->_maskArray[jdx];
          if (pos >= 0 && pos < profileLength)
            candidates[pos] = 0;
        }
      }
    }
  }

  for (size_t idx = 0; idx < this->profileLength; idx++) {
    if (candidates[idx] != 0)
      return true;
  }

  return false;

  //TODO:N's
}


bool CountProfile::checkForSpuriousTransitionDropsWithWindow(uint32_t *maxProfile, unsigned int covEst, double percDrop)
{

  unsigned short kmerSpan = this->translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  unsigned int corrValues[this->profileLength];
  memset(corrValues, 0, sizeof(*corrValues) * this->profileLength);

  // correct values with probability for observing the same sequencing error multiple times
  std::vector<unsigned int> corrWindow;
  double corrFactor = 0.001;

  for(unsigned int idx = 0; idx < kmerSpan/2+1; idx++)
    corrWindow.push_back(this->profile[idx].count);
  corrValues[0] = corrFactor * *(std::max_element(corrWindow.begin(), corrWindow.end())) + 1;

  for(unsigned int idx = 1; idx < this->profileLength; idx++) {
    if (idx + kmerSpan/2 < this->profileLength)
      corrWindow.push_back(this->profile[idx + kmerSpan/2].count);
    if (idx > kmerSpan/2)
      corrWindow.erase(corrWindow.begin());
    unsigned int max = *(std::max_element(corrWindow.begin(), corrWindow.end()));
    corrValues[idx] = (unsigned int)(corrFactor * max+1);

  }

  unsigned int candidates[this->profileLength];
  memset(candidates, 0, sizeof(*candidates) * this->profileLength);

  unsigned int maxProfileLength = this->profileLength + kmerSpan - 1;

  vector<uint32_t> window_back;
  vector<uint32_t> window_front;

  // fill windows for first entry
  window_back.push_back(this->profile[0].count);
  for (size_t idx = 1; idx < WINDOW_SIZE+1; idx++){
    window_front.push_back(this->profile[idx].count);
  }

  uint32_t dropstart_level = UINT_MAX, dropend_level = UINT_MAX, compare_level;
  unsigned int dropstart = this->profileLength, dropend = this->profileLength;

  if(this->profile[0].count < covEst/3.0)
    dropstart = 0;
  for (size_t idx = 1; idx < this->profileLength; idx++){

    uint32_t window_back_level = *(std::min_element(std::begin(window_back), std::end(window_back)));
    uint32_t window_front_level = *(std::min_element(std::begin(window_front), std::end(window_front)));

    if ((double) this->profile[idx].count < percDrop * ((double) window_back_level-(double)corrValues[idx-1]) &&\
        window_back_level - this->profile[idx].count >= SIGNIFICANT_LEVEL_DIFF && \
        this->profile[idx].count < (double)covEst/3.0){

      //drop start
      dropstart = idx;
      dropend = this->profileLength;
      dropstart_level = window_back_level;
    }

    else if (window_front.size() == WINDOW_SIZE && \
            ((percDrop * ((double) window_front_level - (double) corrValues[idx] ) >  (double) window_back_level
              && window_front_level - window_back_level >= SIGNIFICANT_LEVEL_DIFF) || window_front_level > (double)covEst/3.0)){

      //drop end (only considered if we already observed a dropstart or if dropstart=0)
      if (dropstart != this->profileLength) {
        dropend = idx;

        dropend_level = window_front_level;
        compare_level = std::min(dropstart_level, dropend_level);
        // set all positions within drop boundaries as candidates if they are smaller than percDrop * compare_level
        for (size_t d = dropstart; d < dropend; d++) {
          if (this->profile[d].count < percDrop * compare_level) {
            candidates[d] = compare_level;
          }
        }

        // check candidates
        for (size_t d = dropstart; d < dropend; d++) {
          if (maxProfile[d] < percDrop * compare_level && \
             ((d > 0 && maxProfile[d] < percDrop * maxProfile[d - 1]) || maxProfile[d] < percDrop * maxProfile[d + 1])) {
            for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
              int pos = d - translator->_maskArray[jdx];
              if (pos >= dropstart && pos < this->profileLength)
                candidates[pos] = 0;
            }
          }
        }
      }

      // reset values
      dropstart = this->profileLength;
      dropstart_level = UINT_MAX;
      dropend_level = UINT_MAX;
      dropend = this->profileLength;
    }

    // update windows
    if (window_back.size()>= WINDOW_SIZE)
      window_back.erase(window_back.begin());
    window_back.push_back(this->profile[idx].count);

    window_front.erase(window_front.begin());
    if (idx + WINDOW_SIZE < this->profileLength)
      window_front.push_back(this->profile[idx+WINDOW_SIZE].count);
  }

  // check if there is a drop region until the end of the read (dropstart without dropend)
  if(dropstart != 0 && dropstart != this->profileLength){
    for (size_t d = dropstart; d < dropend; d++) {
      if (this->profile[d].count < percDrop * dropstart_level) {
        candidates[d] = dropstart_level;
      }
    }
    for (size_t d = dropstart; d < maxProfileLength; d++) {
      if (maxProfile[d] < percDrop * dropstart_level && \
         ((d > 0 && maxProfile[d] < percDrop * maxProfile[d - 1]) || \
          (d + 1 < maxProfileLength &&  maxProfile[d] < percDrop * maxProfile[d + 1]))) {
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = d - translator->_maskArray[jdx];
          if (pos >= dropstart)
            if(pos < this->profileLength)
              candidates[pos] = 0;
        }
      }
    }
  }

  // check for unexplained drops -> transition drops
  for(unsigned int idx=0; idx< this->profileLength; idx++){
    if (candidates[idx] != 0)
      return true;
  }

  return false;
}