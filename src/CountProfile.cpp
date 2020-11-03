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

  return maxProfile;
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


bool CountProfile::checkForTransitionDrops(unsigned int minCount) {
  unsigned short kmerSpan = translator->getSpan();
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


  char to_mask[maxProfileLen];
  memset(to_mask, 0, maxProfileLen * sizeof(char));

  vector<unsigned int> positions; // to_mask == 2
  vector<unsigned int> candidates; // to_mask >= 1
  int last = -1;
  for (size_t idx = 1; idx < maxProfileLen; idx++) {
    if ((double) (maxProfile[idx]) / maxProfile[idx - 1] < 0.1 && maxProfile[idx] < minCount ) {

      last = 0;
        candidates.push_back(idx);
        to_mask[idx] = 1;


    }
    else if (((double) (maxProfile[idx - 1]) / maxProfile[idx] < 0.1) && maxProfile[idx-1] < minCount) { //} && (last!=-1 || idx>41) ) {
      last = 1;
      if(to_mask[idx - 1] == 0){
      candidates.push_back(idx - 1);
      to_mask[idx - 1] = 1;
      }
    }
  }

  //vector<unsigned int> dropPositions = positions;
  vector<unsigned int> dropPositions = candidates;
  bool checkPoints[profileLength + kmerSpan - 1];
  memset(checkPoints, true, sizeof(bool) * (profileLength + kmerSpan - 1));

  /*
  for (size_t idx = 0; idx < dropPositions.size(); idx++) {
    std::cout << dropPositions[idx] << ",";
  }
  std::cout << std::endl;*/

  for (size_t idx = 0; idx < dropPositions.size(); idx++) {

    for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
      int pos = dropPositions[idx] - translator->_maskArray[jdx];
      if (pos >= 0)
        checkPoints[pos] = false;
    }
  }

  //special mask for N positions
  for (size_t idx = 0; idx < profileLength; idx++) {
    if (!profile[idx].valid)
      checkPoints[idx] = false;
  }


  std::cout << seqinfo->name.c_str() << "\t";
  for (size_t idx = 0; idx < profileLength + kmerSpan - 1; idx++) {
    if (checkPoints[idx] == false)
    std::cout << idx << ",";
  }
  std::cout << std::endl;

  unsigned int dropstart=0,dropend=0;

  for (size_t idx = 1; idx < this->profileLength; idx++) {
    if (//(double) (this->profile[idx].count < minCount) &&
        (double) this->profile[idx].count / this->profile[idx - 1].count < 0.1) {

      if (dropstart == 0 || dropstart == profileLength)
        dropstart = idx;
      dropend = this->profileLength;
    } else if (//(double) (this->profile[idx - 1].count < minCount) &&
               (double) this->profile[idx - 1].count / this->profile[idx].count < 0.1) {
      dropend = idx;
      if ((dropstart !=0 && dropend - dropstart > 0) || dropend-dropstart > 3) {
        unsigned int count = 0;
        for (short jdx = dropstart; jdx < dropend; jdx++) {
          if (checkPoints[jdx] && (double) (this->profile[jdx].count < minCount))
            count++;
        }
        //std::cout << "dropstart: " << dropstart << ", dropend" << dropend<< std::endl;
        //std::cout << "count: " << count << std::endl;
        if (count > 0){// && count <= 41) {
          return true;
        }
        dropstart = this->profileLength;
      }
    }
  }
  if ((double) (this->profile[profileLength-1].count < minCount)) {
    if (dropend - dropstart > 3){//} && dropend - dropstart <= 41) {
      unsigned int count = 0;
      for (short jdx = dropstart; jdx < dropend; jdx++) {
        if (checkPoints[jdx] && (double) (this->profile[jdx].count < minCount))
          count++;
      }
      if (count > 0){//} && count <= 41) {
        return true;
      }
    }
  }

  return false;

}

unsigned int CountProfile::calcMedian() {
  //copy max count values
  uint32_t *max_count = new uint32_t[profileLength];

  for (size_t idx = 0; idx < profileLength; idx++) {
    max_count[idx] = profile[idx].count;
  }

  // sort in ascending order
  sort(max_count, max_count + profileLength);

  auto abundanceEstimation = max_count[(uint32_t) (0.5 * (double) this->profileLength)];

  delete[] max_count;
  return abundanceEstimation;
}
bool CountProfile::checkForTransitionDropsNew(unsigned int minCount) {

  float lowLevelCriterion=0.4;

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  unsigned int median = this->calcMedian();

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

  for (size_t idx = 0; idx < maxProfileLen; idx++)
    std::cout << maxProfile[idx] <<",";
  std::cout << std::endl;
  char to_mask[maxProfileLen];
  memset(to_mask, 0, maxProfileLen * sizeof(char));

  vector<unsigned int> positions; // to_mask == 2
  vector<unsigned int> candidates; // to_mask >= 1
  int last = -1;
  for (size_t idx = 1; idx < maxProfileLen; idx++) {

    if(((double) maxProfile[idx] < lowLevelCriterion*median && (double) maxProfile[idx-1] < lowLevelCriterion*median) ||
      ((double) maxProfile[idx] > lowLevelCriterion*median && (double) maxProfile[idx-1] > lowLevelCriterion*median))
        continue;

    if ((double) (maxProfile[idx]) / maxProfile[idx - 1] < 0.5 && maxProfile[idx] < lowLevelCriterion*median ) {

      candidates.push_back(idx);
      to_mask[idx] = 1;


    }
    else if (((double) (maxProfile[idx - 1]) / maxProfile[idx] < 0.5) && maxProfile[idx-1] < lowLevelCriterion*median ) {
      last = 1;
      if(to_mask[idx - 1] == 0){
        candidates.push_back(idx - 1);
        to_mask[idx - 1] = 1;
      }
    }
  }

  //vector<unsigned int> dropPositions = positions;
  vector<unsigned int> dropPositions = candidates;
  bool checkPoints[profileLength + kmerSpan - 1];
  memset(checkPoints, true, sizeof(bool) * (profileLength + kmerSpan - 1));


  for (size_t idx = 0; idx < dropPositions.size(); idx++) {
    std::cout << dropPositions[idx] << ",";
  }
  std::cout << std::endl;

  for (size_t idx = 0; idx < dropPositions.size(); idx++) {

    for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
      int pos = dropPositions[idx] - translator->_maskArray[jdx];
      if (pos >= 0)
        checkPoints[pos] = false;
    }
  }

  //special mask for N positions
  for (size_t idx = 0; idx < profileLength; idx++) {
    if (!profile[idx].valid)
      checkPoints[idx] = false;
  }


  std::cout << seqinfo->name.c_str() << "\t";
  for (size_t idx = 0; idx < profileLength + kmerSpan - 1; idx++) {
    if (checkPoints[idx] == false)
      std::cout << idx << ",";
  }
  std::cout << std::endl;

  unsigned int dropstart=0,dropend=0;

  for (size_t idx = 1; idx < this->profileLength; idx++) {

    if(((double) this->profile[idx].count < lowLevelCriterion*median && (double) this->profile[idx-1].count < lowLevelCriterion*median) ||
      ((double) this->profile[idx].count > lowLevelCriterion*median && (double) this->profile[idx-1].count > lowLevelCriterion*median))
      continue;

    if (//(double) (this->profile[idx].count < minCount) &&
      (double) this->profile[idx].count / this->profile[idx - 1].count < 0.5) {

      if (dropstart == 0 || dropstart == profileLength)
        dropstart = idx;
      dropend = this->profileLength;
    } else if (//(double) (this->profile[idx - 1].count < minCount) &&
      (double) this->profile[idx - 1].count / this->profile[idx].count < 0.5) {
      dropend = idx;
      if ((dropstart !=0 && dropend - dropstart > 0) || dropend-dropstart > 3) {
        unsigned int count = 0;
        std::cout << "dropstart: " << dropstart << std::endl;
        std::cout << "dropend: " << dropend << std::endl;
        for (short jdx = dropstart; jdx < dropend; jdx++) {
          if (checkPoints[jdx] && (double) (this->profile[jdx].count < lowLevelCriterion*median))
            count++;
        }
        //std::cout << "dropstart: " << dropstart << ", dropend" << dropend<< std::endl;
        //std::cout << "count: " << count << std::endl;
        if (count > 0){// && count <= 41) {
          return true;
        }
        dropstart = this->profileLength;
      }
    }
  }
  if ((double) (this->profile[profileLength-1].count < minCount)) {
    if (dropend - dropstart > 3){//} && dropend - dropstart <= 41) {
      unsigned int count = 0;
      for (short jdx = dropstart; jdx < dropend; jdx++) {
        if (checkPoints[jdx] && (double) (this->profile[jdx].count < lowLevelCriterion*median))
          count++;
      }
      if (count > 0){//} && count <= 41) {
        return true;
      }
    }
  }

  return false;

}


bool CountProfile::checkForTransitionDropsNew2(unsigned int minCount) {

  double lowLevelCriterion = 0.1;
  //double mincount = 10;

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  unsigned int median = this->calcMedian();

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
    std::cout << idx << ": " << maxProfile[idx] << std::endl;
  }
  unsigned int dropstart = 0, dropend = this->profileLength;

  uint8_t candidates[this->profileLength];
  memset(candidates, 0, sizeof(*candidates) * this->profileLength);
  for (size_t idx = 1; idx < this->profileLength; idx++) {

    if (((double) this->profile[idx].count < lowLevelCriterion * median &&
         (double) this->profile[idx - 1].count < lowLevelCriterion * median) ||
        ((double) this->profile[idx].count > lowLevelCriterion * median &&
         (double) this->profile[idx - 1].count > lowLevelCriterion * median))
      continue;

    // dropstart
    if ((double) this->profile[idx].count / this->profile[idx - 1].count < 0.5 &&
        (double) this->profile[idx].count < lowLevelCriterion * median) {
      dropstart = idx;
      dropend = this->profileLength;
    }
      // dropend
    else if ((double) this->profile[idx - 1].count / this->profile[idx].count < 0.5 &&
             (double) this->profile[idx - 1].count < lowLevelCriterion * median) {
      dropend = idx - 1;

      if (dropstart < this->profileLength) {
        std::cout << "dropstart: " << dropstart << std::endl;
        std::cout << "dropend: " << dropend << std::endl;
        for (size_t d = dropstart; d <= dropend; d++) {
          if (this->profile[d].count < lowLevelCriterion * median)
            candidates[d] = 1;
        }

        if (maxProfile[dropend] < lowLevelCriterion * median) {
          std::cout << "Maks position " << dropend << std::endl;
          for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            int pos = dropend - translator->_maskArray[jdx];
            if (pos >= 0)
              candidates[pos] = 0;
          }
        }

        if (maxProfile[dropstart + kmerSpan - 1] < lowLevelCriterion * median) {
          std::cout << "Maks position " << dropstart + kmerSpan - 1 << std::endl;
          for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            int pos = dropstart + kmerSpan - 1 - translator->_maskArray[jdx];
            if (pos >= 0)
              candidates[pos] = 0;
          }
        }
        dropstart = this->profileLength;
      }
      dropend = this->profileLength;
    }
  }


  if (dropstart < this->profileLength) {
    dropend = this->profileLength - 1;
    std::cout << "dropstart: " << dropstart << std::endl;
    std::cout << "dropend: " << dropend << std::endl;
    for (size_t d = dropstart; d <= dropend; d++) {
      if (this->profile[d].count < lowLevelCriterion * median)

        candidates[d] = 1;
    }

    if (maxProfile[dropend] < lowLevelCriterion * median) {
      for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
        int pos = dropend - translator->_maskArray[jdx];
        if (pos >= 0)
          candidates[pos] = 0;
      }
    }

    if (maxProfile[dropstart + kmerSpan - 1] < lowLevelCriterion * median) {
      for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
        int pos = dropstart + kmerSpan - 1 - translator->_maskArray[jdx];
        if (pos >= 0)
          candidates[pos] = 0;
      }
    }
  }


  for (size_t idx = 0; idx < this->profileLength; idx++) {
    if (candidates[idx] != 0) {
      return true;
    }
  }
  //TODO: N's

  return false;

}
//
//bool CountProfile::checkForSpuriousTransitionDrops_old(unsigned int minCount) {
//
//
//  unsigned int median = this->calcMedian();
//  double dropLevelCriterion = 0.01 * median;
//
//  uint8_t candidates[this->profileLength];
//  memset(candidates, 0, sizeof(*candidates) * this->profileLength);
//
//  unsigned int dropstart=0;
//  unsigned int dropend=this->profileLength;
//  for (size_t idx = 1; idx < this->profileLength; idx++) {
//
//    if (((double) this->profile[idx].count < dropLevelCriterion &&
//         (double) this->profile[idx - 1].count < dropLevelCriterion) ||
//        ((double) this->profile[idx].count > dropLevelCriterion &&
//         (double) this->profile[idx - 1].count > dropLevelCriterion))
//      continue;
//
//    // dropstart
//    if ((double) this->profile[idx].count / this->profile[idx - 1].count < 0.5 &&
//        (double) this->profile[idx].count < dropLevelCriterion) {
//      dropstart = idx;
//
//    }
//      // dropend
//    else if ((double) this->profile[idx - 1].count / this->profile[idx].count < 0.5 &&
//             (double) this->profile[idx - 1].count < dropLevelCriterion) {
//      dropend = idx;
//
//      for (size_t d = dropstart; d < dropend; d++) {
//        candidates[d] = 1;
//      }
//      dropstart=this->profileLength;
//      dropend = this->profileLength;
//    }
//  }
//
//  if(dropstart>0) {
//    for (size_t d = dropstart; d < this->profileLength; d++) {
//      candidates[d] = 1;
//    }
//  }
//
//  unsigned short kmerSpan = translator->getSpan();
//  unsigned short kmerWeight = translator->getWeight();
//
//  size_t maxProfileLen = this->profileLength + kmerSpan - 1;
//  uint32_t maxProfile[maxProfileLen];
//  for (size_t idx = 0; idx < maxProfileLen; idx++)
//    maxProfile[idx] = 1;
//
//  for (size_t idx = 0; idx < this->profileLength; idx++) {
//    for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
//      size_t pos = idx + translator->_maskArray[jdx];
//      if (this->profile[idx].valid)
//        maxProfile[pos] = std::max(this->profile[idx].count, maxProfile[pos]);
//    }
//  }
//
//  dropstart=0;
//  dropend=maxProfileLen;
//  //TODO: edge case ?
//  for (size_t idx = 1; idx < maxProfileLen; idx++) {
//    if (((double) maxProfile[idx] < dropLevelCriterion &&
//         (double) maxProfile[idx-1] < dropLevelCriterion) ||
//        ((double) maxProfile[idx] > dropLevelCriterion &&
//         (double) maxProfile[idx-1] > dropLevelCriterion))
//      continue;
//
//    // dropstart
//    if ((double) maxProfile[idx] / maxProfile[idx-1] < 0.5 &&
//        (double) maxProfile[idx] < dropLevelCriterion * median) {
//      dropstart = idx;
//    }
//      // dropend
//    else if ((double) maxProfile[idx-1] / maxProfile[idx] < 0.5 &&
//             (double) maxProfile[idx-1] < dropLevelCriterion * median) {
//      dropend = idx;
//
//      /*for (size_t d = dropstart; d < dropend; d++) {
//        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
//          int pos = d - translator->_maskArray[jdx];
//          if (pos >= 0)
//            candidates[pos] = 0;
//        }
//      }*/
//      for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
//        int pos = dropstart - translator->_maskArray[jdx];
//        if (pos >= 0)
//          candidates[pos] = 0;
//      }
//
//      for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
//        int pos = dropend-1 - translator->_maskArray[jdx];
//        if (pos >= 0)
//          candidates[pos] = 0;
//      }
//
//      dropstart=maxProfileLen;
//      dropend=maxProfileLen;
//    }
//  }
//
//  if(dropstart>0) {
//    /*for (size_t d = dropstart; d < maxProfileLen; d++) {
//      for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
//        int pos = d - translator->_maskArray[jdx];
//        if (pos >= 0)
//          candidates[pos] = 0;
//      }
//    }*/
//    for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
//      int pos = dropstart - translator->_maskArray[jdx];
//      if (pos >= 0)
//        candidates[pos] = 0;
//    }
//
//    for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
//      int pos = dropend-1 - translator->_maskArray[jdx];
//      if (pos >= 0)
//        candidates[pos] = 0;
//    }
//  }
//
//  for (size_t idx = 0; idx < this->profileLength; idx++) {
//    if (candidates[idx] != 0) {
//      return true;
//    }
//  }
//  return false;
//
//  //TODO:N's
//}


char CountProfile::checkForSpuriousTransitionDrops(uint32_t *maxProfile, unsigned int dropLevelCriterion, bool maskOnlyDropEdges) {



  uint8_t candidates[this->profileLength];
  memset(candidates, 0, sizeof(*candidates) * this->profileLength);

  double correctionFactor= 0.01;
  unsigned int dropstart = 0;
  unsigned int dropend = this->profileLength;
  for (size_t idx = 1; idx < this->profileLength; idx++) {

    if (((double) this->profile[idx].count <= dropLevelCriterion &&
         (double) this->profile[idx - 1].count <= dropLevelCriterion) ||
        ((double) this->profile[idx].count >= dropLevelCriterion &&
         (double) this->profile[idx - 1].count >= dropLevelCriterion))
      continue;

    // dropstart
    if ((double) this->profile[idx].count / this->profile[idx - 1].count <= 0.5 &&
        (double) this->profile[idx].count <= dropLevelCriterion) {
      dropstart = idx;

    }
      // dropend
    else if ((double) this->profile[idx - 1].count / this->profile[idx].count <= 0.5 &&
             (double) this->profile[idx - 1].count <= dropLevelCriterion) {

      if (dropstart == this->profileLength)
        continue;

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

  for (size_t idx = 0; idx < this->profileLength; idx++) {
    std::cout << (int) candidates[idx] << std::endl;
  }

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  size_t maxProfileLen = this->profileLength + kmerSpan - 1;

  dropstart = 0;
  dropend = maxProfileLen;

  for (size_t idx = 1; idx < maxProfileLen; idx++) {

    if (((double) maxProfile[idx] <= (dropLevelCriterion + correctionFactor*maxProfile[idx-1]) &&
         (double) maxProfile[idx-1] <= (dropLevelCriterion + correctionFactor*maxProfile[idx])) ||
        ((double) maxProfile[idx] >= (dropLevelCriterion + correctionFactor*maxProfile[idx-1]) &&
         (double) maxProfile[idx-1] >= (dropLevelCriterion + correctionFactor*maxProfile[idx])))
      continue;

    // dropstart
    if ((double) maxProfile[idx] / maxProfile[idx-1] < 0.5 &&
        (double) maxProfile[idx] <= (dropLevelCriterion + correctionFactor*maxProfile[idx-1])) {
      dropstart = idx;
    }
      // dropend
    else if ((double) maxProfile[idx-1] / maxProfile[idx] < 0.5 &&
             (double) maxProfile[idx-1] <= (dropLevelCriterion + correctionFactor*maxProfile[idx])) {


      if (dropstart == maxProfileLen)
        continue;

      dropend = idx;
      std::cout << "maksing Dropstart: " << dropstart << " dropend: " << dropend << std::endl;
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

  std::cout << "after masking" << std::endl;
  for (size_t idx = 0; idx < this->profileLength; idx++) {
    std::cout << (int) candidates[idx] << std::endl;
  }


  for (size_t idx = 0; idx < this->profileLength; idx++) {
    if (candidates[idx] != 0)
      return true;
  }

  return false;

  //TODO:N's
}

unsigned int CountProfile::calc67quantile() {
  //copy max count values
  uint32_t *max_count = new uint32_t[profileLength];

  for (size_t idx = 0; idx < profileLength; idx++) {
    max_count[idx] = profile[idx].count;
  }

  // sort in ascending order
  sort(max_count, max_count + profileLength);

  // 67% quantile
  auto abundanceEstimation = max_count[(uint32_t) (0.67 * (double) this->profileLength)];

  delete[] max_count;
  return abundanceEstimation;
}
