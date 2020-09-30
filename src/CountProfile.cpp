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

std::vector<unsigned int> CountProfile::getDropPointsInMaximzedProfile() {
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
  /*std::cout << "maximized " << this->seqinfo->name << " " << this->seqinfo->comment << std::endl;
  for (size_t idx = 0; idx < maxProfileLen; idx++) {
    std::cout << idx << "\t" << maxProfile[idx] << std::endl;
  }*/

  char to_mask[maxProfileLen];
  memset(to_mask, 0, maxProfileLen * sizeof(char));
  // 0 - not interesting
  // 1 - candidate position for drop or rise
  // 2 - drop or rise position

  vector<unsigned int> positions; // to_mask == 2
  vector<unsigned int> candidates; // to_mask >= 1
  for (size_t idx = 1; idx < maxProfileLen; idx++) {
    if ((double) (maxProfile[idx]) / maxProfile[idx - 1] < 0.1) {
      candidates.push_back(idx);
      to_mask[idx] = 1;
    }
    else if (((double) (maxProfile[idx - 1]) / maxProfile[idx] < 0.1) && to_mask[idx - 1] == 0) {
      candidates.push_back(idx - 1);
      to_mask[idx - 1] = 1;
    }
  }

  unsigned int prev_identified;
  do{

    unsigned int upperLevel, valid;
    prev_identified = positions.size();

    for (unsigned int candidate : candidates) {

      if (to_mask[candidate] == 1) {

        upperLevel = 0;
        valid = 0;
        for (size_t jdx = 1; valid < MIN_UPPER_LEVEL_POSITIONS && jdx < candidate; jdx++) {
          if (to_mask[candidate - jdx] != 2) {
            valid++;

            if ((double) maxProfile[candidate] / maxProfile[candidate - jdx] < 0.1)
              upperLevel++;
            else
              break;
          }
        }
        if (upperLevel >= MIN_UPPER_LEVEL_POSITIONS) {
          to_mask[candidate] = 2;
          positions.push_back(candidate);
          /*if (candidate != maxProfileLen - 1 && to_mask[candidate + 1] == 0) {
            to_mask[candidate + 1] = 1;
            candidates.push_back(candidate + 1);
          }*/
          continue;
        }

        valid = 0;
        upperLevel = 0;
        for (size_t jdx = 1; valid < MIN_UPPER_LEVEL_POSITIONS && jdx < maxProfileLen - candidate; jdx++) {
          if (to_mask[candidate + jdx] != 2) {
            valid++;

            if ((double) maxProfile[candidate] / maxProfile[candidate + jdx] < 0.1)
              upperLevel++;
            else
              break;
          }
        }
        if (upperLevel >= MIN_UPPER_LEVEL_POSITIONS) {
          to_mask[candidate] = 2;
          positions.push_back(candidate);
          /*if (candidate != 0 && to_mask[candidate - 1] == 0) {
            to_mask[candidate - 1] = 1;
            candidates.push_back(candidate - 1);
          }*/
        }
      }
    }
  } while(prev_identified < positions.size());

  return positions;
}





bool CountProfile::getDropPointsSimplified(unsigned int minCount) {
  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  vector<unsigned int> candidates;
  unsigned int dropstart=0,dropend;
  for (size_t idx = 1; idx < this->profileLength; idx++) {
    if ((double) (this->profile[idx].count < minCount) && (double) this->profile[idx].count/this->profile[idx-1].count < 0.1) {

      dropstart=idx;
      dropend = this->profileLength;
    }
    else if ((double) (this->profile[idx-1].count < minCount) && (double) this->profile[idx-1].count/this->profile[idx].count < 0.1) {
      dropend=idx;
      if (dropend-dropstart > 3 ) { //&& dropend-dropstart <= 41 ) {
        return true;
      }
    }
  }
  if ((double) (this->profile[profileLength-1].count < minCount)) {
    if (profileLength - dropstart > 3 ) { //&& profileLength - dropstart <= 41) {
      return true;
    }
  }

  return false;

}

bool CountProfile::getDropPointsSimplified2(unsigned int minCount) {
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
    if ((double) (maxProfile[idx]) / maxProfile[idx - 1] < 0.1) {

      last = 0;
        candidates.push_back(idx);
        to_mask[idx] = 1;


    }
    else if (((double) (maxProfile[idx - 1]) / maxProfile[idx] < 0.1)) { //} && (last!=-1 || idx>41) ) {
      last = 1;
      if(to_mask[idx - 1] == 0){
      candidates.push_back(idx - 1);
      to_mask[idx - 1] = 1;
      }
    }
  }

  if(last == 0 && candidates[candidates.size()-1] >= maxProfileLen-41) {
    candidates.pop_back();
  }

  unsigned int prev_identified;
  do{

    unsigned int upperLevel, valid;
    prev_identified = positions.size();

    for (unsigned int candidate : candidates) {

      if (to_mask[candidate] == 1) {

        upperLevel = 0;
        valid = 0;
        for (size_t jdx = 1; valid < MIN_UPPER_LEVEL_POSITIONS_NEW && jdx < candidate; jdx++) {
          if (to_mask[candidate - jdx] != 2) {
            valid++;

            if ((double) maxProfile[candidate] / maxProfile[candidate - jdx] < 0.1)
              upperLevel++;
            else
              break;
          }
        }
        if (upperLevel >= MIN_UPPER_LEVEL_POSITIONS_NEW) {
          to_mask[candidate] = 2;
          positions.push_back(candidate);
          /*if (candidate != maxProfileLen - 1 && to_mask[candidate + 1] == 0) {
            to_mask[candidate + 1] = 1;
            candidates.push_back(candidate + 1);
          }*/
          continue;
        }

        valid = 0;
        upperLevel = 0;
        for (size_t jdx = 1; valid < MIN_UPPER_LEVEL_POSITIONS_NEW && jdx < maxProfileLen - candidate; jdx++) {
          if (to_mask[candidate + jdx] != 2) {
            valid++;

            if ((double) maxProfile[candidate] / maxProfile[candidate + jdx] < 0.1)
              upperLevel++;
            else
              break;
          }
        }
        if (upperLevel >= MIN_UPPER_LEVEL_POSITIONS_NEW) {
          to_mask[candidate] = 2;
          positions.push_back(candidate);
          /*if (candidate != 0 && to_mask[candidate - 1] == 0) {
            to_mask[candidate - 1] = 1;
            candidates.push_back(candidate - 1);
          }*/
        }
      }
    }
  } while(prev_identified < positions.size());

  vector<unsigned int> dropPositions = positions;
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

  for (size_t idx = 0; idx < profileLength + kmerSpan - 1; idx++) {
    if (checkPoints[idx] == false)
    std::cout << idx << ",";
  }
  std::cout << std::endl;

  unsigned int dropstart=0,dropend=0;

  for (size_t idx = 1; idx < this->profileLength; idx++) {
    if (//(double) (this->profile[idx].count < minCount) &&
        (double) this->profile[idx].count / this->profile[idx - 1].count < 0.1) {

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
        std::cout << "dropstart: " << dropstart << ", dropend" << dropend<< std::endl;
        std::cout << "count: " << count << std::endl;
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


bool CountProfile::getDropPointsSimplified3(unsigned int minCount) {
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
  for (size_t idx = 1; idx < maxProfileLen; idx++) {
    if ((double) (maxProfile[idx]) / maxProfile[idx - 1] < 0.1) {
      candidates.push_back(idx);
      to_mask[idx] = 1;
    }
    else if (((double) (maxProfile[idx - 1]) / maxProfile[idx] < 0.1) && to_mask[idx - 1] == 0) {
      candidates.push_back(idx - 1);
      to_mask[idx - 1] = 1;
    }
  }

  unsigned int prev_identified;
  do{

    unsigned int upperLevel, valid;
    prev_identified = positions.size();

    for (size_t idx = 0; idx < maxProfileLen; idx++) {

      if (to_mask[idx] < 2) {

        upperLevel = 0;
        valid = 0;
        for (size_t jdx = 1; valid < MIN_UPPER_LEVEL_POSITIONS_NEW && jdx < idx; jdx++) {
          if (to_mask[idx - jdx] != 2) {
            valid++;

            if ((double) maxProfile[idx] / maxProfile[idx - jdx] < 0.1)
              upperLevel++;
            else
              break;
          }
        }
        if (upperLevel >= MIN_UPPER_LEVEL_POSITIONS_NEW) {
          to_mask[idx] = 2;
          positions.push_back(idx);
          /*if (candidate != maxProfileLen - 1 && to_mask[candidate + 1] == 0) {
            to_mask[candidate + 1] = 1;
            candidates.push_back(candidate + 1);
          }*/
          continue;
        }

        valid = 0;
        upperLevel = 0;
        for (size_t jdx = 1; valid < MIN_UPPER_LEVEL_POSITIONS_NEW && jdx < maxProfileLen - idx; jdx++) {
          if (to_mask[idx + jdx] != 2) {
            valid++;

            if ((double) maxProfile[idx] / maxProfile[idx + jdx] < 0.1)
              upperLevel++;
            else
              break;
          }
        }
        if (upperLevel >= MIN_UPPER_LEVEL_POSITIONS_NEW) {
          to_mask[idx] = 2;
          positions.push_back(idx);
          /*if (candidate != 0 && to_mask[candidate - 1] == 0) {
            to_mask[candidate - 1] = 1;
            candidates.push_back(candidate - 1);
          }*/
        }
      }
    }
  } while(prev_identified < positions.size());

  vector<unsigned int> dropPositions = positions;
  bool checkPoints[profileLength + kmerSpan - 1];
  memset(checkPoints, true, sizeof(bool) * (profileLength + kmerSpan - 1));

  std::cout << "dropPositions: ";
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

  std::cout << "not checkPoints: ";
  for (size_t idx = 0; idx < profileLength + kmerSpan - 1; idx++) {
    if (checkPoints[idx] == false)
    std::cout << idx << ",";
  }
  std::cout << std::endl;

  unsigned int dropstart=0,dropend;
  for (size_t idx = 1; idx < this->profileLength; idx++) {
    if ((double) (this->profile[idx].count < minCount) &&
        (double) this->profile[idx].count / this->profile[idx - 1].count < 0.1) {

      dropstart = idx;
      dropend = this->profileLength;
    } else if ((double) (this->profile[idx - 1].count < minCount) &&
               (double) this->profile[idx - 1].count / this->profile[idx].count < 0.1) {
      dropend = idx;
      if (dropend - dropstart > 3) {
        unsigned int count = 0;
        for (short jdx = dropstart; jdx < dropend; jdx++) {
          if (checkPoints[jdx])
            count++;
        }
        //std::cout << "dropstart: " << dropstart << ", dropend" << dropend<< std::endl;
        //std::cout << "count: " << count << std::endl;
        if (count > 0 && count <= 41) {
          return true;
        }
      }
    }
  }
  if ((double) (this->profile[profileLength-1].count < minCount)) {
    if (dropend - dropstart > 3 && dropend - dropstart <= 41) {
      unsigned int count = 0;
      for (short jdx = dropstart; jdx < dropend; jdx++) {
        if (checkPoints[jdx])
          count++;
      }
      if (count > 0 && count <= 41) {
        return true;
      }
    }
  }

  return false;

}

bool CountProfile::checkForRiseAndDropPoints(std::vector<unsigned int> dropPositions, unsigned int minCount) {
  unsigned short kmerWeight = translator->getWeight();
  unsigned short kmerSpan = translator->getSpan();
  bool checkPoints[profileLength + kmerSpan - 1];
  memset(checkPoints, true, sizeof(bool) * (profileLength + kmerSpan - 1));

  std::cout << "dropPositions:";
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

  /*std::cout << "#"  << this->seqinfo->name << std::endl;
  for (size_t kdx = 0; kdx < profileLength; kdx++) {
    std::cout << kdx << '\t' << profile[kdx].count << '\t' << (checkPoints[kdx]?"1":"0") << std::endl;
  }*/

  unsigned int upperLevel;
  for (size_t idx = 0; idx < profileLength; idx++) {
    if (checkPoints[idx]) {
      if (profile[idx].count < minCount) {
        upperLevel = 0;
        for (size_t jdx = 1, validCount = 0; validCount < MIN_UPPER_LEVEL_POSITIONS && jdx < idx; jdx++) {
          if (checkPoints[idx - jdx]) {
            validCount++;
            if ((double) profile[idx].count / profile[idx - jdx].count < 0.1)
              upperLevel++;
            else
              break;
          }
        }
        if (upperLevel >= MIN_UPPER_LEVEL_POSITIONS)
          return true;

        upperLevel = 0;
        for (size_t jdx = 1, validCount = 0; validCount < MIN_UPPER_LEVEL_POSITIONS && jdx < profileLength - idx; jdx++) {
          if (checkPoints[idx + jdx]) {
            validCount++;
            if ((double) profile[idx].count / profile[idx + jdx].count < 0.1)
              upperLevel++;
            else
              break;
          }
        }
        if (upperLevel >= MIN_UPPER_LEVEL_POSITIONS)
          return true;

      }
    }
  }
  return false;
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
