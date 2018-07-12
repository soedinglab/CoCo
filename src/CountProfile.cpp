#include "CountProfile.h"


CountProfile::CountProfile(const KmerTranslator *translator,
                           const Lookuptable *lookuptable)
{
  this->translator = translator;
  this->lookuptable = lookuptable;
  this->profile = NULL;
}

CountProfile::CountProfile(const KmerTranslator *translator,
                           const Lookuptable *lookuptable,
                           const SeqType seq, char* readName)
{
  // call constructor
  this->translator = translator;
  this->lookuptable = lookuptable;

  assert(seq.size() > 0);
  this->maxprofileLength = seq.length();
  this->profileLength = seq.length();
  this->profile = new CountProfileEntry[profileLength];
  this->fill(seq, readName);
}

CountProfile::~CountProfile()
{
  delete[] this->profile;
}

void CountProfile::fill(const SeqType seq, const char* readName)
{
  assert(this->lookuptable != NULL);
  assert(this->translator != NULL);

  unsigned short kmerSpan = translator->getSpan();
  unsigned short kmerWeight = translator->getWeight();

  /* coverage s_ni of the spaced k-mer starting at position <i>,
    compute c_i = max{s_{i-j} : b_j=1 and 0 ≤ j <= k-1} */

  size_t seqlen = seq.size();
  assert(seqlen >= kmerSpan);
  this->readName = (char*) readName;
  this->profileLength = seqlen;
  this->populationCoverage = 0;


  // check size of array and create new one in case of too small size
  if(seqlen > maxprofileLength)
  {
    delete[] this->profile;
    this->profile = new CountProfileEntry[seqlen];
    this->maxprofileLength = seqlen;
  }

  // add new read -> store all nucleotids of current read, one per entry in
  // <profiles>
  spacedKmerType kmer = 0, nStore = 0;
  for (size_t idx = 0; idx < seqlen; idx++)
  {
    if (seq[idx] != -1)
    {
      kmer = (kmer << 2) | seq[idx];
      nStore = nStore << 2;

      this->profile[idx].nuc = seq[idx];
    }
    else // found non valid nucleotide
    {
      kmer = kmer << 2;
      nStore = nStore << 2 | 3;
      this->profile[idx].nuc = alphabetSize;
    }

    this->profile[idx].readPos = idx;
    this->profile[idx].count = 0;

    if(idx >= kmerSpan-1)
    {
      // check if spaced k-mer contains 'N'
      if ((nStore & translator->_mask) != 0)
        continue;

      // get abundance for packed k-mer, started at position idx in
      // current read
      kmerType packedKmer = translator->kmer2minPackedKmer(kmer);
      uint32_t count = lookuptable->getCount(packedKmer);
      // update c_i for all previous positions j with b_j=1 and 0 ≤ j ≤ k-1
      // maximizing over abundance value
      for(size_t jdx = 0; jdx < kmerWeight; jdx++)
      {
        size_t pos = idx-(kmerSpan-translator->_maskArray[jdx]-1);
        profile[pos].count = std::max(count, profile[pos].count);
        //printf("profile[%u].count=%u\n", pos, profile[pos].count);
      }
    }
  }
}

size_t CountProfile::calcPopulationCoverage()
{
  //copy max count values
  uint32_t max_count[profileLength];
  for(size_t idx = 0; idx < profileLength; idx++)
  {
    max_count[idx] = profile[idx].count;
  }

  // sort in ascending order
  sort(max_count, max_count+profileLength);

  // 67% quantile
  populationCoverage = max_count[(uint32_t)(0.67*(double)this->profileLength)];

  return populationCoverage;
}
