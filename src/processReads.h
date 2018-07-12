#ifndef PROCESSREADS
#define PROCESSREADS

#include "CountProfile.h"
#include "Lookuptable.h"
#include "KmerTranslator.h"
#include "kseq.h"

KSEQ_INIT(int, read)

FILE* openFileOrDie(const char *fileName, const char * mode)
{

  FILE* file;
  file = fopen(fileName, mode);
  if(file == NULL)
  {
    fprintf(stderr, "ERROR: opening failed for file %s\n", fileName);
    perror(fileName);
    exit(EXIT_FAILURE);
  }
  return file;
}

int processReadFile(const char* readFilename,
                    const char* resultFilename,
                    const Lookuptable &lookuptable,
                    const KmerTranslator &translator)
{
  FILE *readFile = openFileOrDie(readFilename, "r");
  kseq_t *seq = kseq_init(fileno(readFile));

  FILE *resultFile = openFileOrDie(resultFilename, "w");
  CountProfile countprofile(&translator, &lookuptable);
  //TODO LOCAL(countprofile);
  printf("iterate over read\n");
  /* iterate over every single read  */
  while (kseq_read(seq) >= 0)
  {
    const size_t len = seq->seq.l;
    const char* readSeq = seq->seq.s;
    const char* readName = seq->name.s;
    SeqType seqStr; seqStr.reserve(len);

    unsigned int kmerSpan = translator.getSpan();
    if(len < kmerSpan)
    {
      //fprintf(stderr, "WARNING: read %lu is too short!\n", readIdx);
      continue;
    }

    /* sequence to 2bit representation */
    for(unsigned int pos = 0; pos < len; pos++)
    {
      seqStr.push_back((char)res2int[(int)readSeq[pos]]);
    }

    /* fill profile */
    countprofile.fill(seqStr, readName);

    /* estimate population coverage */
    uint32_t popCoverage = countprofile.calcPopulationCoverage();
    fprintf(resultFile, "%s\t%u\n", readName, popCoverage);

  }
  kseq_destroy(seq);
  fclose(readFile);
  fclose(resultFile);
  return EXIT_SUCCESS;
}

#endif // PROCESSREADS

