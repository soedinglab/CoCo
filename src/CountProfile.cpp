// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#include <cstring>

#include "CountProfile.h"
#include "mathsupport.h"
#include "Info.h"
#include "kmer.h"

#define WINDOW_SIZE 4


CountProfile::CountProfile(const KmerTranslator *translator,
                           const LookupTableBase *lookuptable) {
  this->translator = translator;
  this->lookuptable = lookuptable;
  this->profile = NULL;
  this->profile_length = 0;
  this->profile_length_alloc = 0;
}

CountProfile::~CountProfile() {
  delete[] this->profile;
}

void CountProfile::update(){

  unsigned short kmerSpan = translator->getSpan();

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

void CountProfile::fill(SequenceInfo *seqinfo, size_t length) {
  assert(this->lookuptable != NULL);
  assert(this->translator != NULL);

  unsigned short kmerSpan = translator->getSpan();
  //unsigned short kmerWeight = translator->getWeight();

  assert(length >= kmerSpan);
  this->seqinfo = seqinfo;
  this->profile_length = length - kmerSpan + 1;


  this->update();
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
  delete[] maxProfile;
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

bool CountProfile::tryInsertionCorrection(unsigned int insertionStart, unsigned int insertionLen, unsigned int threshold,
                                                                                     uint32_t *neighborhoodTolerance){

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
      assert(res2int[(int) seqinfo->seq[seqPosForFirstKmer]] != -1);
      firstKmer =
        (firstKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForFirstKmer]]);

      unsigned int seqPosForSecondLastKmer = secondLastKmerStart + translator->_mask_array[jdx];
      if (seqPosForSecondLastKmer >= insertionStart)
        seqPosForSecondLastKmer += insertionLen;
      assert(res2int[(int) seqinfo->seq[seqPosForSecondLastKmer]] != -1);
      secondLastKmer =
        (secondLastKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForSecondLastKmer]]);

      unsigned int seqPosForMidKmer = midKmerStart + translator->_mask_array[jdx];
      if (seqPosForMidKmer >= insertionStart)
        seqPosForMidKmer += insertionLen;
      assert(res2int[(int) seqinfo->seq[seqPosForMidKmer]] != -1);
      midKmer =
        (midKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForMidKmer]]);
    }

    // check count change
    unsigned int approved = 0;
    uint32_t c1 = lookuptable->getCount(translator->kmer2minPackedKmer(firstKmer));
    if (c1 > threshold + neighborhoodTolerance[firstKmerStart]) //TODO
      approved++;
    uint32_t c2 = lookuptable->getCount(translator->kmer2minPackedKmer(midKmer));
    if (c2 > threshold + neighborhoodTolerance[midKmerStart]) //TODO
      approved++;
    uint32_t c3 = lookuptable->getCount(translator->kmer2minPackedKmer(secondLastKmer));
    if (c3 > threshold + neighborhoodTolerance[secondLastKmerStart]) //TODO
      approved++;

    if (approved == 3)
      return true;
  }

  return false;
}

int CountProfile::tryDeletionCorrection(unsigned int deletionPos, unsigned int threshold,
                                uint32_t *neighborhoodTolerance){

  // try to fix deletion by inserting a nucleotide at deletionPos

  // deletionPos can neither be the first or last nucleotide of the sequence, because we can not see that
  // so we always find a k-mer which starts before deletionPos and a kmer which contain positions behind deletionPos

  unsigned short kmerSpan = this->translator->getSpan();
  unsigned int firstKmerStart = deletionPos >= kmerSpan? deletionPos-kmerSpan+1:0;
  unsigned int lastKmerStart = deletionPos < profile_length? deletionPos-1:profile_length;
  // attention: profile_length-1 is normally the last valid k-mer, but we will try to insert a nucleotide tehrefore the
  // last valid k-mer starts at profile_length

  packedKmerType firstKmer = 0, lastKmer = 0;
  //
  for (size_t jdx = 0; jdx < translator->weight; jdx++) {

    //TODO: function to build gapped kmer
    unsigned int seqPosForFirstKmer = firstKmerStart + translator->_mask_array[jdx];
    if(seqPosForFirstKmer < deletionPos)
      firstKmer =  (firstKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForFirstKmer]]);
    else if (seqPosForFirstKmer == deletionPos)
      firstKmer =  (firstKmer << 2);
    else
      firstKmer =  (firstKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForFirstKmer-1]]);

    unsigned int seqPosForLastKmer = lastKmerStart + translator->_mask_array[jdx];
    if(seqPosForLastKmer < deletionPos)
      lastKmer =  (lastKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForLastKmer]]);
    else if (seqPosForLastKmer == deletionPos)
      lastKmer =  (lastKmer << 2);
    else
      lastKmer =  (lastKmer << 2) | (3 & res2int[(int) seqinfo->seq[seqPosForLastKmer-1]]);

  }

  int mutationTarget = -1;
  for (unsigned int resMutation = 0; resMutation < alphabetSize; resMutation++) {
    unsigned int improvement = 0;

    //TODO: function to replace nucleotide

    firstKmer &= ~((packedKmerType) 3 << (packedKmerType) (2 * (translator->getWeight() -
                                                                translator->_inverse_mask_array[deletionPos -
                                                                                                firstKmerStart] - 1)));

    firstKmer |= ((packedKmerType) resMutation << (packedKmerType) (2 * (translator->getWeight() -
                                                                         translator->_inverse_mask_array[deletionPos -
                                                                                                         firstKmerStart] -
                                                                         1)));

    lastKmer &= ~((packedKmerType) 3 << (packedKmerType) (2 * (translator->getWeight() -
                                                               translator->_inverse_mask_array[deletionPos -
                                                                                               lastKmerStart] - 1)));

    lastKmer |= ((packedKmerType) resMutation << (packedKmerType) (2 * (translator->getWeight() -
                                                                        translator->_inverse_mask_array[deletionPos -
                                                                                                        lastKmerStart] -
                                                                        1)));

    uint32_t c1 = lookuptable->getCount(translator->kmer2minPackedKmer(firstKmer));
    if (c1 > threshold + neighborhoodTolerance[deletionPos])
      improvement++;
    uint32_t c2 = lookuptable->getCount(translator->kmer2minPackedKmer(lastKmer));
    if (c2 > threshold + neighborhoodTolerance[deletionPos])
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

int CountProfile::doIndelCorrection(uint32_t *maxProfile, unsigned int threshold, double tolerance){
  unsigned short kmerSpan = this->translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  unsigned int maxProfileLength = profile_length + kmerSpan - 1;
  uint32_t neighborhoodTolerance[maxProfileLength];
  memset(neighborhoodTolerance, 0, sizeof(*neighborhoodTolerance) * maxProfileLength);

  // calculate neighborhood tolerance values with probability for observing the same sequencing error multiple times
  std::vector<unsigned int> neighborhood;

  for (unsigned int idx = 0; idx < (unsigned int) kmerSpan / 2 + 1; idx++)
    neighborhood.push_back(maxProfile[idx]);
  neighborhoodTolerance[0] = tolerance * *(std::max_element(neighborhood.begin(), neighborhood.end())) + 1;

  for (unsigned int idx = 1; idx < maxProfileLength; idx++) {
    if (idx + kmerSpan / 2 < maxProfileLength)
      neighborhood.push_back(maxProfile[idx + kmerSpan / 2]);
    if (idx > kmerSpan / 2)
      neighborhood.erase(neighborhood.begin());
    unsigned int max = *(std::max_element(neighborhood.begin(), neighborhood.end()));
    neighborhoodTolerance[idx] = (unsigned int) (tolerance * max + 1);
  }

  //TODO: outsource neighborhood calculation to avoid code duplication and code repetition

  bool changed = false;
  unsigned int dropLen = 0;
  int offset = 0;
  string sequence = seqinfo->seq;
  for (size_t idx = 0; idx <= this->profile_length; idx++) {
    if (idx < profile_length && this->profile[idx].count <= threshold + neighborhoodTolerance[idx]) {
      dropLen++;
    } else if (dropLen > 0) {

      // only check single indel errors
      if (dropLen <= kmerSpan) {

        bool insertion_approved = false, deletion_approved = false;
        unsigned int deletionPos = UINT_MAX, insertionPos = UINT_MAX;
        int resToAdd = -1;

        // set position to correct
        if (idx == profile_length) {// no rise was found, drops with 0 < dropLen <= kmerSpan to end
          deletionPos = idx - dropLen + kmerSpan - 1;
          insertionPos = idx - dropLen + kmerSpan - 1;

        } else if (idx == dropLen || dropLen >= kmerSpan - 1) {
          //  * drops with 0 < dropLen <= kmerSpan from start
          //  * drops kmerSpan-1 <= dropLen <= kmerSpan within the read
          deletionPos = idx ;
          insertionPos = idx - 1;

        }

        // try correction
        if (insertionPos !=UINT_MAX && dropLen > 1)
          insertion_approved = this->tryInsertionCorrection(insertionPos, 1, threshold,
                                                            neighborhoodTolerance);
        if (deletionPos != UINT_MAX && dropLen != kmerSpan) {
          resToAdd = this->tryDeletionCorrection(deletionPos, threshold, neighborhoodTolerance);
          if (resToAdd >= 0)
            deletion_approved = true;
          else if (resToAdd == -2){
            insertion_approved = false;
            // reset insertion_approved because multiple possibilities for deletion correction were found
          }
        }

        // assign correction
        if (insertion_approved && !deletion_approved) {
          // final elimination of insertion error
          sequence.erase(sequence.begin() + offset + insertionPos);
          offset -= 1;
        } else if (deletion_approved && !insertion_approved) {
          // final elimination of deletion error
          sequence.insert(deletionPos + offset, 1, int2res[resToAdd]);
          offset += 1;
        }
        //else if (insertion_approved && deletion_approved){ => both approved => ambigious choice => no correction}
        //else{ no approved => no correction }
        //eliminate error in first/last nucleotide if no correction could be applied
        else if ((idx == 1 || idx == this->profile_length) && dropLen == 1) {

          if (idx == this->profile_length)
            sequence.pop_back();
          else
            sequence.erase(sequence.begin());
          offset -= 1;
        }

      }
      dropLen = 0;
    }
    //else outside drop
  }
  seqinfo->seq = sequence;

  return 0; //TODO
}


int CountProfile::doSubstitutionCorrection(uint32_t *maxProfile, unsigned int covEst, unsigned int threshold, double tolerance, bool dryRun) {

  unsigned short kmerSpan = this->translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();
  unsigned int maxProfileLength = profile_length + kmerSpan - 1;
  uint32_t neighborhoodTolerance[maxProfileLength];
  memset(neighborhoodTolerance, 0, sizeof(*neighborhoodTolerance) * maxProfileLength);

  // calculate neighborhood tolerance values with probability for observing the same sequencing error multiple times
  std::vector<unsigned int> neighborhood;

  for (unsigned int idx = 0; idx < (unsigned int) kmerSpan / 2 + 1; idx++)
    neighborhood.push_back(maxProfile[idx]);
  neighborhoodTolerance[0] = tolerance * *(std::max_element(neighborhood.begin(), neighborhood.end())) + 1;

  for (unsigned int idx = 1; idx < maxProfileLength; idx++) {
    if (idx + kmerSpan / 2 < maxProfileLength)
      neighborhood.push_back(maxProfile[idx + kmerSpan / 2]);
    if (idx > kmerSpan / 2)
      neighborhood.erase(neighborhood.begin());
    unsigned int max = *(std::max_element(neighborhood.begin(), neighborhood.end()));
    neighborhoodTolerance[idx] = (unsigned int) (tolerance * max + 1);

  }

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
    if (maxProfile[idx] <= threshold + neighborhoodTolerance[idx]) {
      errorPositions[foundErrors] = idx;
      for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
        int pos = idx - translator->_mask_array[jdx];
        if (pos >= 0 && pos < (int) profile_length) {
          affectedPostions[pos] = affectedPostions[pos] | ((uint64_t) 1 << (uint64_t) foundErrors);
        }
      }
      foundErrors++;

      if (dryRun) {
        // error candidate positions
        Info(Info::DEBUG) << "error\t" << seqinfo->name.c_str() << "\t" << idx << "\n";
      }
    }
  }
  if (foundErrors == 0) {
    return ERROR_FREE;
  }
  if (dryRun)
    return NONE_CORRECTED;

  //  2. correct sequencing errors
  size_t correctedErrors = 0;

  for (size_t idx = 0; idx < foundErrors; idx++) {

    char current_res = seqinfo->seq[errorPositions[idx]];
    unsigned int start = errorPositions[idx] >= translator->span ? errorPositions[idx] - translator->span + 1 : 0;
    unsigned int end = errorPositions[idx] < profile_length ? errorPositions[idx] : profile_length - 1;
    unsigned int firstUniqueKmerStart = UINT_MAX, lastUniqueKmerStart = 0;
    for (size_t jdx = start; jdx <= end; jdx++) {
      if ((affectedPostions[jdx] > 0) && ((affectedPostions[jdx] & (1 << idx)) == affectedPostions[jdx])) {
        firstUniqueKmerStart = std::min(firstUniqueKmerStart, (unsigned int) jdx);
        lastUniqueKmerStart = std::max(lastUniqueKmerStart, (unsigned int) jdx);
      }
    }
    // first-last-unique-kmer correction strategy
    packedKmerType firstUniqueKmer = 0, lastUniqueKmer = 0;

    // independent error
    if (firstUniqueKmerStart <= lastUniqueKmerStart) {
      for (size_t jdx = 0; jdx < translator->weight; jdx++) {
        firstUniqueKmer =
          (firstUniqueKmer << 2) | (3 & res2int[(int) seqinfo->seq[firstUniqueKmerStart + translator->_mask_array[jdx]]]);
        lastUniqueKmer =
          (lastUniqueKmer << 2) | (3 & res2int[(int) seqinfo->seq[lastUniqueKmerStart + translator->_mask_array[jdx]]]);
      }
    } else {
      continue; //TODO: add strategy for non independent errors
    }
    char mutationTarget = current_res;
    for (unsigned int resMutation = 0; resMutation < alphabetSize; resMutation++) {
      unsigned int improvement = 0;
      if (resMutation == (unsigned int) res2int[(int) current_res])
        continue;

      firstUniqueKmer &= ~((packedKmerType) 3 << (packedKmerType) (2 * (translator->getWeight() -
                                                                        translator->_inverse_mask_array[
                                                                          errorPositions[idx] -
                                                                          firstUniqueKmerStart] - 1)));
      firstUniqueKmer |= ((packedKmerType) resMutation << (packedKmerType) (2 * (translator->getWeight() -
                                                                                 translator->_inverse_mask_array[
                                                                                   errorPositions[idx] -
                                                                                   firstUniqueKmerStart] - 1)));

      lastUniqueKmer &= ~((packedKmerType) 3 << (packedKmerType) (2 * (translator->getWeight() -
                                                                       translator->_inverse_mask_array[
                                                                         errorPositions[idx] -
                                                                         lastUniqueKmerStart] - 1)));
      lastUniqueKmer |= ((packedKmerType) resMutation << (packedKmerType) (2 * (translator->getWeight() -
                                                                                translator->_inverse_mask_array[
                                                                                  errorPositions[idx] -
                                                                                  lastUniqueKmerStart] - 1)));

      uint32_t c1 = lookuptable->getCount(translator->kmer2minPackedKmer(firstUniqueKmer));
      if (c1 > threshold + neighborhoodTolerance[errorPositions[idx]])
        improvement++;
      uint32_t c2 = lookuptable->getCount(translator->kmer2minPackedKmer(lastUniqueKmer));
      if (c2 > threshold + neighborhoodTolerance[errorPositions[idx]])
        improvement++;

      if (improvement == 2) {
        if (mutationTarget == current_res) {
          mutationTarget = int2res[resMutation];
        } else {
          //ambigious mutation choice -> skip correction
          mutationTarget = current_res;
          break;
        }
      }
      // TODO: 1. majority voting? 2. local covEst 3. global covEst (= median)
    }

    if (mutationTarget != current_res) {
      Info(Info::DEBUG) << "corrected\t" << seqinfo->name.c_str() << "\t" << errorPositions[idx]
                        << "\t" << seqinfo->seq[errorPositions[idx]] << "\t" << mutationTarget << "\n";
      seqinfo->seq[errorPositions[idx]] = mutationTarget;
      //TODO: update count table?
      correctedErrors++;
    }
  }


  //Had errors
  if (correctedErrors == 0) {
    return NONE_CORRECTED;
  } else if (correctedErrors < foundErrors) {
    return SOME_CORRECTED;
  } else
    return ALL_CORRECTED;

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

