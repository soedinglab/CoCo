// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include "CountProfile.h"
#include "mathsupport.h"


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

  std::vector<unsigned int> positions;
  for (size_t idx = 1; idx < maxProfileLen; idx++) {
    if ((double) (maxProfile[idx]) / maxProfile[idx - 1] < 0.1)
      positions.push_back(idx);
    else if (((double) (maxProfile[idx - 1]) / maxProfile[idx] < 0.1) &&
             (positions.empty() || positions.back() != idx - 1))
      positions.push_back(idx - 1);
  }

  return positions;
}

bool CountProfile::checkForRiseAndDropPoints(std::vector<unsigned int> dropPositions, unsigned int minCount) {
  unsigned short kmerWeight = translator->getWeight();
  unsigned short kmerSpan = translator->getSpan();
  bool checkPoints[profileLength + kmerSpan - 1];
  memset(checkPoints, true, sizeof(bool) * (profileLength + kmerSpan - 1));


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


  /*
    std::cout << "#"  << this->seqinfo->name << std::endl;
  for (size_t kdx = 0; kdx < profileLength; kdx++) {
    std::cout << kdx << '\t' << profile[kdx].count << '\t' << (checkPoints[kdx]?"1":"0") << std::endl;
  }
   */
  unsigned int validCount = 0;
  size_t idx = 0;
  for (size_t idx = 0; idx < profileLength; idx++) {
    if (checkPoints[idx]) {
      validCount = profile[idx].count;
      idx++;
      break;
    }

  }

  for (; idx < profileLength; idx++) {
    if (checkPoints[idx]) {
      if (((double) profile[idx].count / validCount < 0.1 && profile[idx].count < minCount) ||
          ((double) validCount / profile[idx].count < 0.1 && validCount < minCount))
        return true;

      validCount = profile[idx].count;
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
