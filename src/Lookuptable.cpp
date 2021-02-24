// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <assert.h>
#include <iostream>
#include <string.h>
#include "Lookuptable.h"

#include "Info.h"
#include "Options.h"
#include "util.h"


Lookuptable::Lookuptable(const size_t nbItems, int countMode) {
  // TODO: assert(LOGINDEXSIZE+LOGOFFSETSIZE == 2*ksize);

  indexGridTable = (size_t *) calloc(indexGridTableSize, sizeof(size_t));
  offsetTable = (IndexEntry *) calloc(nbItems, sizeof(IndexEntry));
  maxNumberItems = nbItems;
  numberItems = 0;
  mode = countMode;

  /* set masks */
  _indexmask = (ipow(2, LOGINDEXSIZE) - 1) << LOGOFFSETSIZE;
  _offsetmask = ipow(2, LOGOFFSETSIZE) - 1;
}

Lookuptable::~Lookuptable() {
  free(this->indexGridTable);
  free(this->offsetTable);
}

inline size_t Lookuptable::getGridPosition(packedKmerType kmer) const {
  return (kmer & _indexmask) >> LOGOFFSETSIZE;
}

inline size_t Lookuptable::getOffset(packedKmerType kmer) const {
  return (kmer & _offsetmask);
}

void Lookuptable::assignKmertoGrid(packedKmerType kmer) {
  size_t gridPosition = getGridPosition(kmer);
  indexGridTable[gridPosition]++;
}


void Lookuptable::setupIndexGridTable() {
  size_t prevElementLength = indexGridTable[0];
  indexGridTable[0] = 0;
  for (size_t i = 0; i < indexGridTableSize - 1; i++) {
    const size_t currElementLength = indexGridTable[i + 1];
    indexGridTable[i + 1] = indexGridTable[i] + prevElementLength;
    prevElementLength = currElementLength;
  }
}

size_t Lookuptable::addElement(packedKmerType kmer, unsigned int count) {
  size_t prevWritingPosition = 0;
  const size_t gridPosition = getGridPosition(kmer);
  const size_t writingPosition = indexGridTable[gridPosition];
  const size_t offset = getOffset(kmer);

  /* if kmer already in lookuptable, increase only count field */
  if (gridPosition > 0)
    prevWritingPosition = indexGridTable[gridPosition - 1];
  if (writingPosition >= this->maxNumberItems) {
    Info(Info::ERROR) << "ERROR: Lookuptable addElement overflows. Current writing position is "
                      << writingPosition << "\n";
    EXIT(EXIT_FAILURE);
  }
  for (size_t pos = prevWritingPosition; pos < writingPosition; pos++) {
    if (offsetTable[pos].count != 0 && offsetTable[pos].indexOffset == offset) {
      if (mode == Options::COUNT_MODE_SUM)
        offsetTable[pos].count += count;
      else if (mode == Options::COUNT_MODE_MAX)
        offsetTable[pos].count = std::max(offsetTable[pos].count,count);
      return pos;
    }
  }

  /* add new kmer in lookuptable */
  offsetTable[writingPosition].indexOffset = offset;
  offsetTable[writingPosition].count = count;
  indexGridTable[gridPosition]++;
  numberItems++;

  return writingPosition;
}

void Lookuptable::finalSetupTables(size_t countThreshold) {

  size_t prev = 0,
    readpos = 0,
    writepos = 0;

  for (size_t idx = 0; idx < indexGridTableSize; idx++) {
    for (; readpos < indexGridTable[idx]; readpos++) {

      if (offsetTable[readpos].count > countThreshold) {
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
  offsetTable = (IndexEntry *) realloc(offsetTable, maxNumberItems * sizeof(IndexEntry));
}

inline std::pair<size_t, size_t> Lookuptable::getIndexGridRange(packedKmerType kmer) const {
  size_t gridPosition = getGridPosition(kmer);
  size_t indexGridSize = indexGridTable[gridPosition + 1] -
                         indexGridTable[gridPosition];
  return std::make_pair(indexGridTable[gridPosition], indexGridSize);
}

unsigned int Lookuptable::getCount(const packedKmerType kmer) const {

  const std::pair<size_t, size_t> grid = getIndexGridRange(kmer);
  const size_t kmerOffset = getOffset(kmer);
  size_t pos = grid.first;
  size_t endPos = pos + grid.second;
  //use != instead of < to save sort step
  for (; pos < endPos && (offsetTable[pos].indexOffset != kmerOffset
                          || offsetTable[pos].count == 0); pos++) {
  }

  return (pos != endPos) ? offsetTable[pos].count : 0;
}

void Lookuptable::iterateOverAll(FILE *fp) const{
  size_t readpos = 0;
  for (size_t idx = 0; idx < indexGridTableSize-1; idx++) {
    char *prefix = packedKmer2String(idx, LOGINDEXSIZE/2);
    for (; readpos < indexGridTable[idx+1]; readpos++) {

      char *suffix = packedKmer2String(offsetTable[readpos].indexOffset, LOGOFFSETSIZE/2);
      fprintf(fp, "%s%s\t%d\n", prefix, suffix, offsetTable[readpos].count);
      free(suffix);
    }
    free(prefix);
  }
}