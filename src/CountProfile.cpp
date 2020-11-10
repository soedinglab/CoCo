// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include "CountProfile.h"
#include "mathsupport.h"

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
      std::cout << "candidate drop start: " << dropstart << " dropend " << dropend << std::endl;
      dropstart=this->profileLength;
      dropend = this->profileLength;
    }
  }

  if(dropstart>0 && dropstart < this->profileLength) {
    std::cout << "candidate dropstart: " << dropstart << " dropend: " << dropend << std::endl;
    for (size_t d = dropstart; d < this->profileLength; d++) {
      if (profile[d].count <= dropLevelCriterion)
        candidates[d] = 1;
    }
  }


  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  size_t maxProfileLen = this->profileLength + kmerSpan - 1;
std::cout << "maxprofile: " << std::endl;
  for(int i=0; i<maxProfileLen;i++)
{
std::cout << i << ": " << maxProfile[i] << std::endl;
}
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


      std::cout << "masking dropstart: " << dropstart << " dropend: " << dropend << std::endl;
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
std::cout << "masking dropstart: " << dropstart << " dropend: " << dropend << std::endl;
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
      {
std::cout << "Remaingin candidate: " << idx << std::endl;
return true;}
  }

  return false;

  //TODO:N's
}

char CountProfile::checkForSpuriousTransitionDropsSupported(uint32_t *maxProfile, unsigned int dropLevelCriterion, bool maskOnlyDropEdges) {



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


      unsigned int window_width=10;
      vector<uint32_t> window;
      vector<uint32_t> window2;

      for (size_t d=0; d<window_width; d++){
        if(idx>d)
          window.push_back(idx-d);
        if(idx+d<this->profileLength)
          window2.push_back(idx+d);
      }
      unsigned int window_level=*(std::min_element(std::begin(window), std::end(window)));
      unsigned int window2_level=*(std::min_element(std::begin(window2), std::end(window2)));

      if (this->profile[idx].count / (double) window_level <= 0.5)
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
      std::cout << "candidate drop start: " << dropstart << " dropend " << dropend << std::endl;
      dropstart=this->profileLength;
      dropend = this->profileLength;
    }
  }

  if(dropstart>0 && dropstart < this->profileLength) {
    std::cout << "candidate dropstart: " << dropstart << " dropend: " << dropend << std::endl;
    for (size_t d = dropstart; d < this->profileLength; d++) {
      if (profile[d].count <= dropLevelCriterion)
        candidates[d] = 1;
    }
  }


  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  size_t maxProfileLen = this->profileLength + kmerSpan - 1;
  std::cout << "maxprofile: " << std::endl;
  for(int i=0; i<maxProfileLen;i++)
  {
    std::cout << i << ": " << maxProfile[i] << std::endl;
  }
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


      std::cout << "masking dropstart: " << dropstart << " dropend: " << dropend << std::endl;
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
    std::cout << "masking dropstart: " << dropstart << " dropend: " << dropend << std::endl;
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
    {
      std::cout << "Remaingin candidate: " << idx << std::endl;
      return true;}
  }

  return false;

  //TODO:N's
}


char CountProfile::checkForSpuriousTransitionDropsWindow(uint32_t *maxProfile, unsigned int dropLevelCriterion, double perc, bool maskOnlyDropEdges) {

  std::cout << "checkForSpuriousTransitionDropsWindow" << std::endl;
  unsigned int candidates[this->profileLength];
  memset(candidates, 0, sizeof(*candidates) * this->profileLength);

 //Todo: try median and mark on the fly instead of minimum and backward marking?
  unsigned int window_width=10;
  vector<uint32_t> window;
  vector<uint32_t> window2;

  window.push_back(this->profile[0].count);

  for (size_t idx = 1; idx < window_width+1; idx++){
    window2.push_back(this->profile[idx].count);
  }

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  size_t maxProfileLen = this->profileLength + kmerSpan - 1;

  uint32_t window_level,window2_level,dropstart_level,dropend_level,compare_level;
  unsigned int dropstart = this->profileLength, dropend;
  for (size_t idx = 1; idx < this->profileLength; idx++){

    window_level=*(std::min_element(std::begin(window), std::end(window)));
    window2_level=*(std::min_element(std::begin(window2), std::end(window2)));
    std::cout << idx << "\t" << window_level << "\t" << window2_level << std::endl;
    if ((double) this->profile[idx].count < 0.5* (double)window_level && this->profile[idx].count < perc *dropLevelCriterion){ //drop start

      std::cout << "found possible dropstart" << idx << std::endl;
      if(dropstart == this->profileLength) {
        std::cout << "set new dropstart" << idx << std::endl;
        dropstart = idx;
        dropend = this->profileLength;
        // std::cout << "dropstart: " << dropstart << std::endl;
        dropstart_level = window_level;
      }
    }
   // else if ( (double) this->profile[idx].count > 2* (double)window_level && this->profile[idx].count == window2_level){ //drop end
    else if ( (double) 0.5*window2_level >  (double)window_level && window2_level > 10){ //drop end
      std::cout << "found possible dropend" << idx << std::endl;
      if (dropstart != this->profileLength) {
        //std::cout << "dropend: " << dropend << " without start"<< std::endl;
        //continue;
        //}

        dropend = idx;
        //std::cout << "dropend: " << dropend << std::endl;
        dropend_level = this->profile[idx].count;

        compare_level = std::min(dropstart_level, dropend_level);
        for (size_t d = dropstart; d < dropend; d++) {
          if (this->profile[d].count < 0.5 * compare_level) {
            candidates[d] = dropstart_level;

            //candidates[d] = compare_level;
            /*if (maxProfile[d] < 0.5 * compare_level) {
              //to_mask
              for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
                int pos = d - translator->_maskArray[jdx];
                if (pos >= 0 && pos < profileLength)
                  candidates[pos] = 0;
              }

            }*/
          }
        }

        for (size_t d = dropstart; d < dropend; d++) {
          if (candidates[d] > 0){
            for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
              int pos = d + translator->_maskArray[jdx];
              if((pos < dropend || d + kmerSpan > this->profileLength) && candidates[d]*0.5 > maxProfile[pos]){
                candidates[d] = 0;
                break;
              }
            }
          }
        }

          unsigned int tmp=0;
          for (size_t d = dropstart; d < dropend; d++) {
            if (candidates[d] >0)
              tmp++;
          }

        std::cout << this->seqinfo->name << "\t" << dropstart << "\t" << dropend << "\t"
                  << dropstart_level << "\t" << dropend_level << "\t" << tmp << "\t"
                  << this->profileLength << "\t" << dropLevelCriterion
                  << std::endl;

        dropstart = this->profileLength;
      } /*else
      {
        std::cout << this->seqinfo->name << "\t" << 0 << "\t" << idx << "\t"
                  << 0 << "\t" << this->profile[idx].count << "\t" << idx << "\t"
                  << this->profileLength
                  << std::endl;
      }*/

    }
    if (window.size()>= window_width)
      window.erase(window.begin());

    window.push_back(this->profile[idx].count);

    window2.erase(window2.begin());
    if (idx+window_width < this->profileLength)
      window2.push_back(this->profile[idx+window_width].count);
  }

  /*for(size_t d= this->profileLength; d<maxProfileLen; d++){
    if (maxProfile[d] < 0.5 * compare_level) {
      //to_mask
      for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
        int pos = d - translator->_maskArray[jdx];
        if (pos >= 0 && pos < profileLength)
          candidates[pos] = 0;
      }
  }*/
  /*for(size_t d=0; d<this->profileLength; d++)
    {
      if (candidates[d] > 0){
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = d + translator->_maskArray[jdx];
          if(candidates[d]*0.5 > maxProfile[pos]){
            candidates[d] = 0;
            break;
          }
        }
      }
    }*/
  /*for(size_t d=0; d<maxProfileLen; d++){
      if (maxProfile[d] < 0.5*candidates[d]){
    }
  }*/

/*
  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  size_t maxProfileLen = this->profileLength + kmerSpan - 1;

  window.clear();
  window.push_back(maxProfile[0]);
  dropstart = maxProfileLen;
  for (size_t idx = 1; idx < maxProfileLen; idx++){

    level=*(std::min_element(std::begin(window), std::end(window)));
    if (maxProfile[idx] < 0.5*level){ //drop start
      dropstart = idx;
      dropstart_level = maxProfile[idx];
    }
    else if (maxProfile[idx]> 0.5*level){ //drop end

      if (dropstart == maxProfileLen)
        continue;

      dropend = idx;
      dropend_level = maxProfile[idx];

      compare_level = std::min(dropstart_level,dropend_level);
      for (size_t d = dropstart; d < dropend; d++) {
        if (maxProfile[d])
      }

      dropstart = maxProfileLen;

    }
    if (window.size()>= window_width)
      window.erase(window.begin());
    window.push_back(maxProfile[idx]);
  }



  for(int i=0; i<maxProfileLen;i++)
  {
    std::cout << i << ": " << maxProfile[i] << std::endl;
  }
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


      std::cout << "masking dropstart: " << dropstart << " dropend: " << dropend << std::endl;
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
    std::cout << "masking dropstart: " << dropstart << " dropend: " << dropend << std::endl;
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
    {
      std::cout << "Remaingin candidate: " << idx << std::endl;
      return true;}
  }

  return false;

  //TODO:N's
  */
}

