#include "CountProfile.h"
#include "mathsupport.h"


CountProfile::CountProfile(const KmerTranslator *translator,
                           const Lookuptable *lookuptable)
{
  this->translator = translator;
  this->lookuptable = lookuptable;
  this->profile = NULL;
}

CountProfile::CountProfile(const KmerTranslator *translator,
                           const Lookuptable *lookuptable,
                           const SeqType seq, char* seqName)
{
  // call constructor
  this->translator = translator;
  this->lookuptable = lookuptable;

  assert(seq.size() > 0);
  this->maxprofileLength = seq.length();
  this->profileLength = seq.length();
  this->profile = new CountProfileEntry[profileLength];
  this->fill(seq, seqName);
}

CountProfile::~CountProfile()
{
  delete[] this->profile;
}

void CountProfile::fill(const SeqType seq, const char* seqName)
{
  assert(this->lookuptable != NULL);
  assert(this->translator != NULL);

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  /* coverage s_ni of the spaced k-mer starting at position <i>,
    compute c_i = max{s_{i-j} : b_j=1 and 0 ≤ j <= k-1} */

  size_t seqlen = seq.size();
  assert(seqlen >= kmerSpan);
  this->seqName = (char*) seqName;
  this->profileLength = seqlen-kmerSpan+1;
  this->abundanceEstimation = 0;


  // check size of array and create new one in case of too small size
  if(this->profileLength > this->maxprofileLength)
  {
    delete[] this->profile;
    this->profile = new CountProfileEntry[this->profileLength];
    this->maxprofileLength = this->profileLength;
  }

  // add new sequence -> store all nucleotids of current sequence, one per entry in
  // <profiles>
  spacedKmerType kmer = 0, nStore = 0;
  for (size_t idx = 0; idx < seqlen; idx++)
  {
    if (seq[idx] != -1)
    {
      kmer = (kmer << 2) | seq[idx];
      nStore = nStore << 2;

      //this->profile[idx].nuc = seq[idx];
    }
    else // found non valid nucleotide
    {
      kmer = kmer << 2;
      nStore = nStore << 2 | 3;
      //this->profile[idx].nuc = alphabetSize;
    }

    /*this->profile[idx].seqPos = idx;
    this->profile[idx].count = 0;*/

    if(idx >= kmerSpan-1)
    {
      // check if spaced k-mer contains 'N'
      if ((nStore & translator->_mask) != 0)
        continue;

      // get abundance for packed k-mer, started at position idx in
      // current sequnece
      kmerType packedKmer = translator->kmer2minPackedKmer(kmer);
      uint32_t count = lookuptable->getCount(packedKmer);
      profile[idx-(kmerSpan-1)].count = count;
      // update c_i for all previous positions j with b_j=1 and 0 ≤ j ≤ k-1
      // maximizing over abundance value
      /*for(size_t jdx = 0; jdx < kmerWeight; jdx++)
      {
        size_t pos = idx-(kmerSpan-translator->_maskArray[jdx]-1);
        profile[pos].count = std::max(count, profile[pos].count);
        
      }*/
    }
  }
}

size_t CountProfile::calc67quantile()
{
  //copy max count values
  uint32_t *max_count = new uint32_t[profileLength];
  
  for(size_t idx = 0; idx < profileLength; idx++)
  {
    max_count[idx] = profile[idx].count;
  }

  // sort in ascending order
  sort(max_count, max_count+profileLength);

  // 67% quantile
  abundanceEstimation = max_count[(uint32_t)(0.67*(double)this->profileLength)];

  delete[] max_count;
  return abundanceEstimation;
}

inline uint32_t getMedian(multiset<uint32_t> &multiSet, size_t size)
{
    double a = *next(multiSet.begin(), size / 2 - 1);
    double b = *next(multiSet.begin(), size / 2);
    if(size & 1)
        return(b);

    return((a + b) * 0.5);
}

inline uint32_t getMax(multiset<uint32_t> &multiSet, size_t size)
{
    return(*next(multiSet.begin(), size-1));
}

std::vector<unsigned int> CountProfile::getDropPointsInMaximzedProfile()
{
    unsigned short kmerSpan = translator->getSpan();
    unsigned short kmerWeight = translator->getWeight();

    size_t maxProfileLen = this->profileLength + kmerSpan - 1;
    uint32_t maxProfile[maxProfileLen];
    memset(maxProfile,0,sizeof(uint32_t)*maxProfileLen);

    for (size_t idx = 0; idx < this->profileLength; idx++) {
        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            size_t pos = idx + translator->_maskArray[jdx];
            maxProfile[pos] = std::max(this->profile[idx].count, maxProfile[pos]);
        }
    }

    std::vector<unsigned int>  positions;
    for (size_t idx = 1; idx < maxProfileLen; idx++) {
      if ((double)(maxProfile[idx])/maxProfile[idx-1] < 0.1)
          positions.push_back(idx);
      else if (((double)(maxProfile[idx-1])/maxProfile[idx] < 0.1) && (positions.empty() || positions.back() != idx-1))
          positions.push_back(idx-1);
    }

    return positions;
}

bool CountProfile::checkForRiseAndDropPoints (std::vector<unsigned int> dropPositions, unsigned int windowsize)
{
    unsigned short kmerWeight = translator->getWeight();
    unsigned short kmerSpan = translator->getSpan();
    bool checkPoints[profileLength+kmerSpan-1];
    memset(checkPoints, true, sizeof(bool)*(profileLength+kmerSpan-1));


    for (size_t idx=0; idx < dropPositions.size(); idx ++) {

        for (size_t jdx = 0; jdx < kmerWeight; jdx++) {
            int pos = dropPositions[idx] - translator->_maskArray[jdx];
            if(pos >=0)
                checkPoints[pos] =  false;
        }
    }

    /*std:cerr << "toIgnore: ";
    for (size_t idx=0; idx < profileLength+kmerSpan-1; idx ++) {
        if (!checkPoints[idx])
          std::cerr <<  idx << ",";
    }

    std::cerr<<endl;*/
    unsigned int validCount = 0;
    size_t idx = 0;
    for (size_t idx = 0; idx < profileLength; idx ++) {
        if (checkPoints[idx]) {
            validCount=profile[idx].count;
            idx++;
            break;
        }

    }

    for (;idx < profileLength; idx ++) {
        /*if ((checkPoints[idx] && checkPoints[idx-1]) && (((double) (profile[idx].count)/profile[idx-1].count < 0.1) ||
                (double) (profile[idx-1].count)/profile[idx].count < 0.1)) {
            std::cerr << seqName << " break at " << idx << " " << profile[idx].count << " " << profile[idx-1].count <<std::endl;
            return true;
        }*/
        if (checkPoints[idx]){
            if (((double) profile[idx].count/validCount  < 0.1) || ((double) validCount/profile[idx].count < 0.1))
              return true;

            validCount = profile[idx].count;
        }

    }
    return false;
}

//bool CountProfile::checkForRiseAndDropPoints(unsigned int windowsize)
//{
//    unsigned long checkPoint = ULONG_LONG_MAX;
//    multiset <uint32_t> window;
//    double prevMedian, median;
//
//    //windowsize = translator->getSpan();
//    for(size_t idx = 0; idx < windowsize; idx++)
//    {
//        window.insert(profile[idx].count);
//    }
//    prevMedian = getMax(window, windowsize);
//    //prevMedian = getMedian(window, windowsize);
//    // compare median of sliding window */
//    unsigned short m = 0.2;
//    for(size_t idx = windowsize; idx < profileLength; idx++)
//    {
//        window.erase(window.find(profile[idx-windowsize].count));
//        window.insert(profile[idx].count);
//        median = getMax(window, windowsize);
//        //median = getMedian(window, windowsize);
//        //std::cout << "median: " << median << std::endl;
//
//        /* identify drop or rise */
//        if (median/(prevMedian-m*sqrt(prevMedian)) < 0.1 || prevMedian/(median-m*sqrt(median)) < 0.1)
//        {
////std::cout << "idx: " << idx << std::endl;
////std::cout << "prevMedian: " << prevMedian << std::endl;
////std::cout << "median: " << median << std::endl;
//            return true;
//        }
//        prevMedian = median;
//    }
//
//    return false;
//}

void CountProfile::showProfile(FILE *fp) const
{
  for (size_t idx=0; idx < profileLength; idx++)
  {
    //fprintf(fp,"%u\t%u\n", profile[idx].seqPos, profile[idx].count);
      fprintf(fp,"%u\t%u\n", idx, profile[idx].count);
  }
}
