// Written by Annika Jochheim <annika.jochheim@mpinat.mpg.de>
#include <cstring>

#include "CountProfile.h"
#include "mathsupport.h"
#include "Info.h"
#include "kmer.h"

#define WINDOW_SIZE 5

static void calcNeighborhoodTolerance(uint32_t *neighborhoodTolerance, uint32_t *maxProfile, unsigned int maxProfileLength,
                                      unsigned short kmerSpan, double threshold, unsigned int pseudocount, unsigned int lowerbound)
{
  memset(neighborhoodTolerance, 0, sizeof(*neighborhoodTolerance) * maxProfileLength);

  // calculate neighborhood pseudocount values with probability for observing the same sequencing error multiple times
  std::vector<unsigned int> neighborhood;

  for (unsigned int idx = 0; idx <= (unsigned int) kmerSpan; idx++)
    neighborhood.push_back(maxProfile[idx]);

  for (unsigned int idx = 0; idx < maxProfileLength; idx++) {
    if (idx > kmerSpan / 2 && idx + kmerSpan / 2 < maxProfileLength) {
      neighborhood.push_back(maxProfile[idx + kmerSpan / 2]);
      neighborhood.erase(neighborhood.begin());
    }

    unsigned int max = *(std::max_element(neighborhood.begin(), neighborhood.end()));

    if (max >= lowerbound) {
      // round(threshold*neighborhood)+pseudocount
      neighborhoodTolerance[idx] = lround(threshold * max) + pseudocount;
    } else
      neighborhoodTolerance[idx] = 0;
  }
}



CountProfile::CountProfile(const KmerTranslator *translator, LookupTableBase *lookuptable) {
  this->translator = translator;
  this->lookuptable = lookuptable;
  this->profile = NULL;
  this->profile_length = 0;
  this->profile_length_alloc = 0;
  this->seqinfo = NULL;
  this->ambigCorr=0;
}

CountProfile::~CountProfile() {
  delete[] this->profile;
}

CountProfileEntry * CountProfile::getProfilecopy() {


  CountProfileEntry *copy = new CountProfileEntry[this->profile_length];
  memcpy(copy, this->profile, this->profile_length*(sizeof(*copy)));
  return copy;
}

void CountProfile::update(bool updateLookuptable){

  if (updateLookuptable){
    for (unsigned int i=0; i < this->profile_length; i++)
      if (this->profile[i].valid)
        this->lookuptable->decreaseCount(this->profile[i].kmer);
  }

  unsigned short kmerSpan = translator->getSpan();
  profile_length = seqinfo->seq.length()- kmerSpan + 1;

  // check size of array and create new one in case of too small size
  if (profile_length > profile_length_alloc) {

    delete[] this->profile;
    this->profile = new CountProfileEntry[profile_length];
    this->profile_length_alloc = profile_length;
  }

  spacedKmerType kmer = 0, nStore = 0;
  for (size_t idx = 0; idx < seqinfo->seq.length(); idx++) {
    if ((char) res2int[(int) seqinfo->seq[idx]] != -1) {
      kmer = (kmer << 2) | (char) res2int[(int) seqinfo->seq[idx]];
      nStore = nStore << 1;
    } else // found non valid nucleotide
    {
      kmer = kmer << 2;
      nStore = nStore << 1 | 1;
    }

    if ((int)idx >= kmerSpan - 1) {
      // check if k-mer contains 'N' (or any other invalid character)
      if ((nStore & translator->_spaced_mask) != 0) {
        profile[idx - (kmerSpan - 1)].count = 0;
        profile[idx - (kmerSpan - 1)].valid = 0;
        profile[idx - (kmerSpan - 1)].kmer = 0;
        continue;
      }

      packedKmerType packedKmer = translator->kmer2minPackedKmer(kmer);
      uint32_t count = lookuptable->getCount(packedKmer);
      profile[idx - (kmerSpan - 1)].count = count;
      profile[idx - (kmerSpan - 1)].valid = 1;
      profile[idx - (kmerSpan - 1)].kmer = packedKmer;
      if (updateLookuptable)
        this->lookuptable->increaseCount(packedKmer);

    }
  }

  //TODO: continue until maxprofillength to reset count and valid flag? not necessary

}

void CountProfile::setSeqInfo(SequenceInfo *seqinfo) {
  assert(this->lookuptable != NULL);
  assert(this->translator != NULL);

  this->seqinfo = seqinfo;
  this->profile_length = 0;
}

void CountProfile::fill(SequenceInfo *seqinfo) {
  assert(this->lookuptable != NULL);
  assert(this->translator != NULL);

  unsigned short kmerSpan = translator->getSpan();
  //unsigned short kmerWeight = translator->getWeight();

  assert(seqinfo->seq.length() >= kmerSpan);
  this->seqinfo = seqinfo;
  this->profile_length = seqinfo->seq.length() - kmerSpan + 1;


  this->update(false);
}

/*** show functions ***/

void CountProfile::showProfile(FILE *fp) const {
  for (size_t idx = 0; idx < profile_length; idx++) {
    fprintf(fp, "%lu\t%u\n", idx, profile[idx].count);
  }
}

void CountProfile::showMaximzedProfile(FILE *fp) const {

  uint32_t *maxProfile = maximize();
  unsigned short kmerSpan = translator->getSpan();

  for (size_t idx = 0; idx < profile_length + kmerSpan - 1; idx++) {
    fprintf(fp, "%lu\t%u\n", idx, maxProfile[idx]);
  }
  delete[] maxProfile;
}

/*** basic profile operations ***/

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

/*** advanced profile operations ***/

int CountProfile::doSubstitutionCorrection(uint32_t *maxProfile, double threshold, unsigned int pseudocount,
                                unsigned int lowerbound, bool needMultipleKmers, bool updateLookup, unsigned int *correctedSubstitutions) {

  unsigned short kmerSpan = this->translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  unsigned int maxProfileLength = profile_length + kmerSpan - 1;
  uint32_t neighborhoodTolerance[maxProfileLength];
  calcNeighborhoodTolerance(neighborhoodTolerance, maxProfile, maxProfileLength, kmerSpan, threshold, pseudocount, lowerbound);

  // 1. find sequencing errors
  unsigned int foundErrors = 0;
  unsigned int errorPositions[64];
  uint64_t affectedPostions[profile_length];
  memset(affectedPostions, 0, sizeof(*affectedPostions) * profile_length);
  memset(errorPositions, 0, sizeof(*errorPositions) * 64);

  for (size_t idx = 0; idx < maxProfileLength; idx++) {
    if (foundErrors == 63) {
      return TOO_MANY_ERRORS;
    }
    if (maxProfile[idx] <= neighborhoodTolerance[idx]) {
      // error candidate positions
      errorPositions[foundErrors] = idx;
      for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
        int pos = idx - translator->_mask_array[jdx];
        if (pos >= 0 && pos < (int) profile_length) {
          affectedPostions[pos] = affectedPostions[pos] | ((uint64_t) 1 << (uint64_t) foundErrors);
        }
      }
      foundErrors++;
    }
  }
  if (foundErrors == 0) {
    return ERROR_FREE;
  }

  //  2. correct sequencing errors
  size_t correctedErrors = 0;
  for (size_t idx = 0; idx < foundErrors; idx++) {

    unsigned int start = errorPositions[idx] >= translator->span ? errorPositions[idx] - translator->span + 1 : 0;
    unsigned int end = errorPositions[idx] < profile_length ? errorPositions[idx] : profile_length - 1;
    unsigned int firstUniqueKmerStart = UINT_MAX, lastUniqueKmerStart = 0;
    for (size_t jdx = start; jdx <= end; jdx++) {
      if ((affectedPostions[jdx] > 0) && (affectedPostions[jdx] == ((uint64_t)1 << idx))) {
        firstUniqueKmerStart = std::min(firstUniqueKmerStart, (unsigned int) jdx);
        lastUniqueKmerStart = std::max(lastUniqueKmerStart, (unsigned int) jdx);
      }
    }

    if(needMultipleKmers && (firstUniqueKmerStart == lastUniqueKmerStart))
      continue; // Skip this error if needMultipleKmers is set, but we only have one kmer

    int mutationTarget = firstLastUniqueKmerCorrectionStrategy(errorPositions[idx], firstUniqueKmerStart,
                                                               lastUniqueKmerStart, neighborhoodTolerance);
    // do correction
    if (mutationTarget >= 0) {
      Info(Info::CDEBUG) << "corrected\t" << seqinfo->name.c_str() << "\t" << errorPositions[idx]
                        << "\t" << seqinfo->seq[errorPositions[idx]] << "\t" << int2res[mutationTarget] << "\n";
      seqinfo->seq[errorPositions[idx]] = int2res[mutationTarget];
      correctedErrors++;
      if(!seqinfo->qual.empty()) {
        seqinfo->qual[errorPositions[idx]] = getAvgQual(seqinfo->qual, errorPositions[idx]);
      }
    }
  }

  (*correctedSubstitutions) += correctedErrors;

  //Had errors
  if (correctedErrors == 0) {
    return NONE_CORRECTED;
  } else if (correctedErrors < foundErrors) {
    this->update(updateLookup);
    return SOME_CORRECTED;
  } else {
    this->update(updateLookup);
    return ALL_CORRECTED;
  }
}

// do indel or edge substitution correction
bool CountProfile::doIndelCorrection(uint32_t *maxProfile, double threshold, unsigned int pseudocount, unsigned int lowerbound, bool trySubstitution,
                                     bool updateLookup, unsigned int *correctedSubstitutions, unsigned int *correctedInsertions, unsigned int *correctedDeletions){
  unsigned short kmerSpan = this->translator->getSpan();
  unsigned int maxProfileLength = profile_length + kmerSpan - 1;
  uint32_t neighborhoodTolerance[maxProfileLength];
  calcNeighborhoodTolerance(neighborhoodTolerance, maxProfile, maxProfileLength, kmerSpan, threshold, pseudocount, lowerbound);

  bool changed = false;
  unsigned int dropLen = 0;
  int offset = 0;
  string sequence = seqinfo->seq;
  string qual = seqinfo->qual;
  // search for drops in count profile
  for (size_t idx = 0; idx <= this->profile_length; idx++) {
    if (idx < profile_length && this->profile[idx].count <= neighborhoodTolerance[idx]) {
      dropLen++;
    } else if (dropLen > 0) { // dropEnd

      // only check single indel or edge substitution errors
      if (dropLen <= kmerSpan) {

        bool insertion_approved = false, deletion_approved = false;
        unsigned int deletionPos = UINT_MAX, insertionPos = UINT_MAX, substitutionPos = UINT_MAX;

        // set positions to correct
        if (idx == profile_length) {
          // reach read end but no rise was found, drops with 0 < dropLen <= kmerSpan at the end
          deletionPos = idx - dropLen + kmerSpan - 1;
          insertionPos = idx - dropLen + kmerSpan - 1;
          substitutionPos = idx - dropLen + kmerSpan - 1;

        } else if (idx == dropLen) {
          // drops with 0 < dropLen <= kmerSpan from start
          deletionPos = idx ;
          insertionPos = idx - 1;
          substitutionPos = idx-1;

        } else if (dropLen >= (unsigned int) kmerSpan - 1) {
          // drops kmerSpan-1 <= dropLen <= kmerSpan within the read
          if (dropLen == (unsigned int) kmerSpan -1)
            deletionPos = idx;
          insertionPos = idx - 1;
        }


        // edge case for substitution error
        if(trySubstitution && (substitutionPos!=UINT_MAX)) {
          int resToSub = this->edgeSubstitutionCorrection(substitutionPos, neighborhoodTolerance);
          if (resToSub >= 0) {
            sequence[substitutionPos + offset] = int2res[resToSub];
            if(!qual.empty()) {
              qual[substitutionPos + offset] = getAvgQual(qual, substitutionPos + offset);
            }
            changed=true;
            (*correctedSubstitutions) +=1;
            continue;
          }
        }
        // try insertion correction
        int resToAdd = -1;
        if (insertionPos != UINT_MAX && dropLen > 1)
          insertion_approved = this->tryInsertionCorrection(insertionPos, 1, neighborhoodTolerance);
        // try deletion correction
        if (deletionPos != UINT_MAX && dropLen != kmerSpan) {
          resToAdd = this->tryDeletionCorrection(deletionPos, neighborhoodTolerance);
          if (resToAdd >= 0)
            deletion_approved = true;
          else if (resToAdd == -2){
            insertion_approved = false;
            // reset insertion_approved because multiple possibilities for deletion correction were found
          }
        }

        // assign indel correction
        if (insertion_approved && !deletion_approved) {
          // final elimination of insertion error
          sequence.erase(sequence.begin() + offset + insertionPos);
          if(!qual.empty()) {
            qual.erase(qual.begin() + offset + insertionPos);
          }
          changed=true;
          offset -= 1;
          (*correctedInsertions) += 1;
        } else if (deletion_approved && !insertion_approved && resToAdd >= 0) {
          // final elimination of deletion error
          sequence.insert(deletionPos + offset, 1, int2res[resToAdd]);
          if(!qual.empty()) {
            qual.insert(deletionPos + offset, 1, 33);
            qual[deletionPos + offset] = getAvgQual(qual, deletionPos + offset);
          }
          changed=true;
          offset += 1;
          (*correctedDeletions) +=1;

        }
        //else if (insertion_approved && deletion_approved){ => both approved => ambigious choice => no correction}
        //else{ no approved => no correction }
        //eliminate error in first/last nucleotide if no correction could be applied
        /*else if ((idx == 1 || idx == this->profile_length) && dropLen == 1) {

          if (idx == this->profile_length)
            sequence.pop_back();
          else
            sequence.erase(sequence.begin());
          offset -= 1;
        }*/

      }
      dropLen = 0;
    }
    //else outside drop
  }
  seqinfo->seq = sequence;
  seqinfo->qual = qual;

  if(changed)
    this->update(updateLookup);

  return changed;
}

bool CountProfile::doTrimming(uint32_t *maxProfile, double threshold, unsigned int pseudocount, unsigned int lowerbound,
                              unsigned int maxTrimLen, bool updateLookup, unsigned int *trimmedCounter){

  bool changed = false;
  unsigned short kmerSpan = this->translator->getSpan();
  unsigned int maxProfileLength = profile_length + kmerSpan - 1;
  uint32_t neighborhoodTolerance[maxProfileLength];
  calcNeighborhoodTolerance(neighborhoodTolerance, maxProfile, maxProfileLength, kmerSpan, threshold, pseudocount,  lowerbound);

  unsigned int dropLen = 0;
  string sequence = seqinfo->seq;
  string qual = seqinfo->qual;


  // look for non corrected drops at start or end of the profile
  for (unsigned int idx = 0; idx < this->profile_length && this->profile[idx].count <= neighborhoodTolerance[idx]; idx++) {
      dropLen++;
  }
  if(dropLen > 0 && dropLen <= maxTrimLen) {
    sequence.erase(sequence.begin(), sequence.begin() + dropLen);
    if (!qual.empty())
      qual.erase(qual.begin(), qual.begin() + dropLen);
    changed = true;
    (*trimmedCounter) += dropLen;
  }

  dropLen = 0;
  for (unsigned int idx = this->profile_length; idx > 0 && this->profile[idx-1].count <= neighborhoodTolerance[idx-1]; idx--) {
    dropLen++;
  }
  if(dropLen > 0 && dropLen <= maxTrimLen ) {
    sequence.erase(sequence.end() - dropLen, sequence.end());
    if (!qual.empty())
      qual.erase(qual.end() - dropLen, qual.end());
    changed = true;
    (*trimmedCounter) += dropLen;
  }

  seqinfo->seq = sequence;
  seqinfo->qual = qual;

  if(changed)
    this->update(updateLookup);

  return changed;
}

int CountProfile::firstLastUniqueKmerCorrectionStrategy(unsigned int errorPos, unsigned int firstUniqueKmerStart,
                                                        unsigned int lastUniqueKmerStart, const uint32_t *neighborhoodTolerance) {
  int mutationTarget = -1;
  char current_res = seqinfo->seq[errorPos];
  packedKmerType firstUniqueKmer = 0, lastUniqueKmer = 0;
  if (firstUniqueKmerStart <= lastUniqueKmerStart) {
    for (size_t jdx = 0; jdx < translator->weight; jdx++) {
      firstUniqueKmer =
              (firstUniqueKmer << 2) | (3 & res2int[(int) seqinfo->seq[firstUniqueKmerStart + translator->_mask_array[jdx]]]);
      lastUniqueKmer =
              (lastUniqueKmer << 2) | (3 & res2int[(int) seqinfo->seq[lastUniqueKmerStart + translator->_mask_array[jdx]]]);
    }
  } else {
    return -1;
  }

  for (unsigned int resMutation = 0; resMutation < alphabetSize; resMutation++) {
    unsigned int improvement = 0;
    if (resMutation == (unsigned int) res2int[(int) current_res])
      continue;

    firstUniqueKmer &= ~((packedKmerType) 3 << (packedKmerType) (2 * (translator->getWeight() -
                                                                      translator->_inverse_mask_array[errorPos - firstUniqueKmerStart] - 1)));
    firstUniqueKmer |= ((packedKmerType) resMutation << (packedKmerType) (2 * (translator->getWeight() -
                                                                               translator->_inverse_mask_array[errorPos - firstUniqueKmerStart] - 1)));

    lastUniqueKmer &= ~((packedKmerType) 3 << (packedKmerType) (2 * (translator->getWeight() -
                                                                     translator->_inverse_mask_array[errorPos - lastUniqueKmerStart] - 1)));
    lastUniqueKmer |= ((packedKmerType) resMutation << (packedKmerType) (2 * (translator->getWeight() -
                                                                              translator->_inverse_mask_array[errorPos - lastUniqueKmerStart] - 1)));

    uint32_t c1 = lookuptable->getCount(translator->kmer2minPackedKmer(firstUniqueKmer));
    if (c1 > neighborhoodTolerance[errorPos])
      improvement++;
    uint32_t c2 = lookuptable->getCount(translator->kmer2minPackedKmer(lastUniqueKmer));
    if (c2 > neighborhoodTolerance[errorPos])
      improvement++;

    if (improvement == 2) {
      if (mutationTarget == -1) {
        mutationTarget = resMutation;
      } else {
        //ambigious mutation choice -> skip correction
        mutationTarget = -2;
        if (firstUniqueKmerStart != lastUniqueKmerStart)
          this->ambigCorr+=1;
        break;
      }
    }
    // TODO: 1. majority voting? 2. local covEst 3. global covEst (= median)
  }

  return mutationTarget;
}

int CountProfile::edgeSubstitutionCorrection(unsigned int substitutionStart, uint32_t *neighborhoodTolerance){


  unsigned int start = substitutionStart >= translator->span ? substitutionStart - translator->span + 1 : 0;
  unsigned int end = substitutionStart < profile_length ? substitutionStart : profile_length - 1;
  unsigned int firstUniqueKmerStart = UINT_MAX, lastUniqueKmerStart = 0;
  //
  for (size_t jdx = start; jdx <= end; jdx++) {
    if (translator->_inverse_mask_array[substitutionStart-jdx] != UCHAR_MAX) {
      //kmer starting at jdx has a coding position at substitutionPos
      firstUniqueKmerStart = std::min(firstUniqueKmerStart, (unsigned int) jdx);
      lastUniqueKmerStart = std::max(lastUniqueKmerStart, (unsigned int) jdx);
    }
  }

  int mutationTarget = firstLastUniqueKmerCorrectionStrategy(substitutionStart, firstUniqueKmerStart, lastUniqueKmerStart,
                                                             neighborhoodTolerance);
  return mutationTarget;

}

bool CountProfile::tryInsertionCorrection(unsigned int insertionStart, unsigned int insertionLen, const uint32_t *neighborhoodTolerance){

  // try to fix insertion by deleting nucleotides from insertionStart to insertionStart+1

  unsigned short kmerSpan = this->translator->getSpan();

  /* choose k-mers for evaluation
     firstKmer = leftmost kmer containing the whole insertion
     secondLastKmer = second rightmost kmer containing the whole insertion
     midKmer = kmer containing whole insertion while counterbalance the number of nucleotides before and after the insertion
  */
  unsigned int firstKmerStart = insertionStart >= kmerSpan? insertionStart-kmerSpan+1:0,
          secondLastKmerStart = insertionStart < profile_length? insertionStart-1: profile_length-insertionLen - 1 ,
          midKmerStart = (firstKmerStart+secondLastKmerStart)/2;

  // make sure we have at least one nucleotide before and after the insert range
  if(firstKmerStart < insertionStart && secondLastKmerStart + kmerSpan >= insertionStart+insertionLen) {

    packedKmerType firstKmer = 0, midKmer = 0, secondLastKmer = 0;
    // remove whole insertion from all three kmers
    for (size_t jdx = 0; jdx < translator->weight; jdx++) {

      unsigned int seqPosForFirstKmer = firstKmerStart + translator->_mask_array[jdx];
      if (seqPosForFirstKmer >= insertionStart)
        seqPosForFirstKmer += insertionLen;
      if(res2int[(int) seqinfo->seq[seqPosForFirstKmer]] == -1)
        return false;
      firstKmer =
              (firstKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForFirstKmer]]);

      unsigned int seqPosForSecondLastKmer = secondLastKmerStart + translator->_mask_array[jdx];
      if (seqPosForSecondLastKmer >= insertionStart)
        seqPosForSecondLastKmer += insertionLen;
      if(res2int[(int) seqinfo->seq[seqPosForSecondLastKmer]] == -1)
        return false;
      secondLastKmer =
              (secondLastKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForSecondLastKmer]]);

      unsigned int seqPosForMidKmer = midKmerStart + translator->_mask_array[jdx];
      if (seqPosForMidKmer >= insertionStart)
        seqPosForMidKmer += insertionLen;
      if(res2int[(int) seqinfo->seq[seqPosForMidKmer]] == -1)
        return false;
      midKmer =
              (midKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForMidKmer]]);
    }

    // check count change
    unsigned int approved = 0;
    uint32_t c1 = lookuptable->getCount(translator->kmer2minPackedKmer(firstKmer));
    if (c1 > neighborhoodTolerance[firstKmerStart])
      approved++;
    uint32_t c2 = lookuptable->getCount(translator->kmer2minPackedKmer(midKmer));
    if (c2 > neighborhoodTolerance[midKmerStart])
      approved++;
    uint32_t c3 = lookuptable->getCount(translator->kmer2minPackedKmer(secondLastKmer));
    if (c3 > neighborhoodTolerance[secondLastKmerStart])
      approved++;

    if (approved == 3)
      return true;
  }

  return false;
}

int CountProfile::tryDeletionCorrection(unsigned int deletionPos, const uint32_t *neighborhoodTolerance){

  // try to fix deletion by inserting a nucleotide at deletionPos

  // deletionPos can neither be the first or last nucleotide of the sequence, because we can not see that
  // so we always find a k-mer which starts before deletionPos and a kmer which contain positions behind deletionPos

  unsigned short kmerSpan = this->translator->getSpan();
  unsigned int firstKmerStart = deletionPos >= kmerSpan? deletionPos-kmerSpan+1:0;
  unsigned int lastKmerStart = deletionPos < profile_length? deletionPos-1:profile_length;
  // attention: profile_length-1 is normally the last valid k-mer, but we will try to insert a nucleotide therefore the
  // last valid k-mer starts at profile_length

  packedKmerType firstKmer = 0, lastKmer = 0;
  //
  for (size_t jdx = 0; jdx < translator->weight; jdx++) {

    unsigned int seqPosForFirstKmer = firstKmerStart + translator->_mask_array[jdx];
    if(seqPosForFirstKmer < deletionPos) {
      if(res2int[(int) seqinfo->seq[seqPosForFirstKmer]] == -1)
        return -1;
      firstKmer = (firstKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForFirstKmer]]);
    }else if (seqPosForFirstKmer == deletionPos)
      firstKmer =  (firstKmer << 2);
    else {
      if(res2int[(int) seqinfo->seq[seqPosForFirstKmer-1]] == -1)
        return -1;
      firstKmer = (firstKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForFirstKmer - 1]]);
    }

    unsigned int seqPosForLastKmer = lastKmerStart + translator->_mask_array[jdx];
    if(seqPosForLastKmer < deletionPos) {
      if (res2int[(int) seqinfo->seq[seqPosForLastKmer]] == -1)
        return -1;
      lastKmer = (lastKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForLastKmer]]);
    }else if (seqPosForLastKmer == deletionPos)
      lastKmer =  (lastKmer << 2);
    else {
      if (res2int[(int) seqinfo->seq[seqPosForLastKmer-1]] == -1)
        return -1;
      lastKmer = (lastKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForLastKmer - 1]]);
    }
  }

  int mutationTarget = -1;
  for (unsigned int resMutation = 0; resMutation < alphabetSize; resMutation++) {
    unsigned int improvement = 0;

    firstKmer &= ~((packedKmerType) 3 << (packedKmerType) (2 * (translator->getWeight() -
                                                                translator->_inverse_mask_array[deletionPos - firstKmerStart] - 1)));

    firstKmer |= ((packedKmerType) resMutation << (packedKmerType) (2 * (translator->getWeight() -
                                                                         translator->_inverse_mask_array[deletionPos - firstKmerStart] - 1)));

    lastKmer &= ~((packedKmerType) 3 << (packedKmerType) (2 * (translator->getWeight() -
                                                               translator->_inverse_mask_array[deletionPos - lastKmerStart] - 1)));

    lastKmer |= ((packedKmerType) resMutation << (packedKmerType) (2 * (translator->getWeight() -
                                                                        translator->_inverse_mask_array[deletionPos - lastKmerStart] -  1)));

    uint32_t c1 = lookuptable->getCount(translator->kmer2minPackedKmer(firstKmer));
    if (c1 > neighborhoodTolerance[deletionPos])
      improvement++;
    uint32_t c2 = lookuptable->getCount(translator->kmer2minPackedKmer(lastKmer));
    if (c2 > neighborhoodTolerance[deletionPos])
      improvement++;

    if (improvement == 2) {
      if (mutationTarget == -1)
        mutationTarget = resMutation;
      else {
        //ambigious mutation choice -> skip correction
        mutationTarget = -2;
        break;
      }
    }
  }
  return mutationTarget;
}

bool CountProfile::checkForSpuriousTransitionDropsWithWindowNew(double threshold)
{
  unsigned short longestBlock = translator->getLongestBlock();
  unsigned short windowSize = std::max(longestBlock+1, WINDOW_SIZE);
  unsigned int span = translator->getSpan();
  std::vector<unsigned int> precidingWindow; precidingWindow.reserve(windowSize);
  std::vector<unsigned int> successiveWindow; successiveWindow.reserve(windowSize);
  int dropstart = 0, dropend;

  for (size_t idx = 0; idx < windowSize; idx++) {
    precidingWindow.push_back(this->profile[idx+windowSize].count);
    successiveWindow.push_back(this->profile[idx].count);
  }
  for (size_t idx = windowSize-1; idx < profile_length-windowSize; idx++){

    successiveWindow.erase(successiveWindow.begin());
    successiveWindow.push_back(this->profile[idx].count);
    precidingWindow.erase(precidingWindow.begin());
    precidingWindow.push_back(this->profile[idx+windowSize].count);

    unsigned int precidingWindow_max = *(std::max_element(precidingWindow.begin(), precidingWindow.end()));
    unsigned int precidingWindow_min = *(std::min_element(precidingWindow.begin(), precidingWindow.end()));

    unsigned int successiveWindow_max = *(std::max_element(successiveWindow.begin(), successiveWindow.end()));
    unsigned int successiveWindow_min = *(std::min_element(successiveWindow.begin(), successiveWindow.end()));

    if(precidingWindow_max == 0 || precidingWindow_min == 0 || successiveWindow_max == 0 || successiveWindow_min == 0)
      continue;

    if((double)precidingWindow_max/(double)successiveWindow_min < threshold) {
      // drop start
      dropstart = idx + 1;
    }
    if((double)successiveWindow_max/(double)precidingWindow_min < threshold) {
      // drop end
      dropend = idx;
      if (dropstart >=0 && dropend-dropstart+1 < span) {
        //std::cout << "found spurious drop" << std::endl;
        return true;
      }
      dropstart = -1;
    }
  }
  if (dropstart >=0 && profile_length-dropstart < span) {
    return true;
  }

  return false;
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
  double corrFactor = 0.001; // TODO:

  for(unsigned int idx = 0; idx < (unsigned int)kmerSpan/2+1; idx++)
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
  //TODO: improve masking as soon coco correction is finished
  //TODO: add parameter to also mask for indels?

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
              if (pos >= (int)dropstart && pos < (int)profile_length)
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
          if (pos >= (int)dropstart)
            if(pos < (int)profile_length)
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
          if (pos >= 0 && pos < (int)profile_length)
            candidates[pos] = 0;
        }

        if (dropend > dropstart + 1) {
          for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            int pos = dropend - 1 - translator->_mask_array[jdx];
            if (pos >= 0 && pos < (int)profile_length)
              candidates[pos] = 0;
          }
        }
      }
      else{
        for (size_t d = dropstart; d < dropend; d++) {
          for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            int pos = d - translator->_mask_array[jdx];
            if (pos >= 0 && pos < (int)profile_length)
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
        if (pos >= 0 && pos < (int)profile_length)
          candidates[pos] = 0;
      }

      if (dropend > dropstart + 1) {
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = dropend - 1 - translator->_mask_array[jdx];
          if (pos >= 0 && pos < (int)profile_length)
            candidates[pos] = 0;
        }
      }
    }
    else {
      for (size_t d = dropstart; d < maxProfileLen; d++) {
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
          int pos = d - translator->_mask_array[jdx];
          if (pos >= 0 && pos < (int)profile_length)
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

