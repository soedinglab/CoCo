#include <assert.h>
#include <iostream>
#include <string.h>
#include "Lookuptable.h"
#include "util.h"


Lookuptable::Lookuptable(const size_t nbItems, float corrFactor)
{
  // TODO: assert(LOGINDEXSIZE+LOGOFFSETSIZE == 2*ksize);

  indexGridTable = (size_t *) calloc(indexGridTableSize, sizeof(size_t));
  offsetTable = (IndexEntry *) calloc(nbItems, sizeof(IndexEntry));
  maxNumberItems = nbItems;
  numberItems = 0;
  this->corrFactor = corrFactor;

  /* set masks */
  _indexmask = (ipow(2,LOGINDEXSIZE)-1) << LOGOFFSETSIZE;
  _offsetmask = ipow(2,LOGOFFSETSIZE)-1;
}

Lookuptable::~Lookuptable()
{
  free(this->indexGridTable);
  free(this->offsetTable);
}

inline size_t Lookuptable::getGridPosition(kmerType kmer) const
{
  return (kmer & _indexmask) >> LOGINDEXSIZE;
}

inline size_t Lookuptable::getOffset(kmerType kmer) const
{
  return (kmer & _offsetmask);
}


void Lookuptable::assignKmertoGrid(kmerType kmer)
{
  size_t gridPosition = getGridPosition(kmer);
  indexGridTable[gridPosition]++;
}


void Lookuptable::setupIndexGridTable()
{
  size_t prevElementLength = indexGridTable[0];
  indexGridTable[0] = 0;
  for(size_t i = 0; i < indexGridTableSize-1; i++)
  {
    const size_t currElementLength = indexGridTable[i + 1];
    indexGridTable[i + 1] = indexGridTable[i] + prevElementLength;
    prevElementLength = currElementLength;
  }
}

size_t Lookuptable::addElement(kmerType kmer, unsigned int count)
{
  size_t prevWritingPosition=0;
  const size_t gridPosition = getGridPosition(kmer);
  const size_t writingPosition = indexGridTable[gridPosition];
  const size_t offset = getOffset(kmer);
  if(writingPosition >= this->maxNumberItems)
  {
    std::cerr << "Lookuptable addElement overflows. Current write position is "
              << writingPosition << std::endl;
    EXIT(EXIT_FAILURE);
  }

  /* if kmer already in lookuptable, increase only count field */
  if (gridPosition >0)
    prevWritingPosition = indexGridTable[gridPosition-1];
  for (size_t pos = prevWritingPosition; pos < writingPosition; pos++)
  {
    if (offsetTable[pos].indexOffset == offset)
    {
      offsetTable[pos].count += count;
      return pos;
    }
  }

  /* add new kmer in lookuptable */
  offsetTable[writingPosition].indexOffset = offset;
  offsetTable[writingPosition].count  = count;
  indexGridTable[gridPosition]++;
  numberItems++;

  return writingPosition;
}

void Lookuptable::finalSetupTables(size_t countThreeshold)
{

  size_t prev = 0,
         readpos = 0,
         writepos = 0;

  for(size_t idx = 0; idx < indexGridTableSize; idx++)
  {
    for(;readpos < indexGridTable[idx]; readpos++)
    {

      if (offsetTable[readpos].count > countThreeshold)
      {
        if (readpos != writepos)
          offsetTable[writepos] = offsetTable[readpos];

        writepos++;
      }
    }
    indexGridTable[idx] = prev;
    prev = writepos;
  }

  assert(writepos <= numberItems);
  numberItems = writepos;
  maxNumberItems = numberItems;
  offsetTable = (IndexEntry *) realloc(offsetTable, maxNumberItems*sizeof(IndexEntry));
}

inline std::pair<size_t, size_t> Lookuptable::getIndexGridRange(kmerType kmer) const
{
    size_t gridPosition = getGridPosition(kmer);
    size_t indexGridSize = indexGridTable[gridPosition + 1] -
                           indexGridTable[gridPosition];
    return std::make_pair(indexGridTable[gridPosition], indexGridSize);
}

unsigned int Lookuptable::getCount (const kmerType kmer) const
{

  const std::pair<size_t, size_t> grid = getIndexGridRange(kmer);
  const size_t kmerOffset = getOffset(kmer);

  size_t pos = grid.first;
  size_t endPos = pos + grid.second;
  //use != instead of < to save sort step
  for(; pos < endPos && (offsetTable[pos].indexOffset != kmerOffset
                         || offsetTable[pos].count == 0); pos++ ){
  }

  return (pos != endPos) ? offsetTable[pos].count : 0;
}

float Lookuptable::getCorrFactor() const
{
  return corrFactor;
}
