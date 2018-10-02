#include "processSequences.h"
#include "CountProfile.h"
#include "Lookuptable.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "util.h"

KSEQ_INIT(int, read)

FILE* openFileOrDie(std::string fileName, const char * mode)
{
  return openFileOrDie(fileName.c_str(), mode);
}

FILE* openFileOrDie(const char *fileName, const char * mode)
{

  FILE* file;
  file = fopen(fileName, mode);
  if(file == NULL)
  {
    fprintf(stderr, "ERROR: opening failed for file %s\n", fileName);
    perror(fileName);
    EXIT(EXIT_FAILURE);
  }
  return file;
}



int processSeqFile(const char* seqFilename,
                   const char* resultFilename,
                   const Lookuptable &lookuptable,
                   const KmerTranslator &translator,
                   int (*processCountProfile)(CountProfile &, FILE*))
{
  FILE *seqFile = openFileOrDie(seqFilename, "r");
  kseq_t *seq = kseq_init(fileno(seqFile));

  FILE *resultFile = openFileOrDie(resultFilename, "w");
  CountProfile countprofile(&translator, &lookuptable);

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
      fprintf(stderr, "WARNING: read %s is too short!\n", readName);
      continue;
    }

    /* sequence to 2bit representation */
    for(unsigned int pos = 0; pos < len; pos++)
    {
      seqStr.push_back((char)res2int[(int)readSeq[pos]]);
    }

    /* fill profile */
    countprofile.fill(seqStr, readName);

    processCountProfile(countprofile, resultFile);
    //TODO: function pointer for what to do with countProfile

  }
  kseq_destroy(seq);
  fclose(seqFile);
  fclose(resultFile);
  return EXIT_SUCCESS;
}

int processSeqFile(string seqFilename,
                   string resultFilename,
                   const Lookuptable &lookuptable,
                   const KmerTranslator &translator,
                   int (*processCountProfile)(CountProfile &, FILE*))
{
  return processSeqFile(seqFilename.c_str(), resultFilename.c_str(),
                        lookuptable, translator,processCountProfile);
}

void process_sampleList(vector<string> *sampleList, std::string resultFileName,
                        Lookuptable *lookuptable,\
                        KmerTranslator *translator,
                        int (*processCountProfile)(CountProfile &, FILE*))
{
  for (vector<string>::iterator sampleIt = sampleList->begin() ; sampleIt != sampleList->end(); ++sampleIt)
  {

    int retval = processSeqFile(sampleIt->c_str(),
                                 resultFileName.c_str(),
                                 *lookuptable,
                                 *translator,
                                 processCountProfile);
    if (retval != EXIT_SUCCESS)
    {
      std::cerr << "ERROR processing read file " << *sampleIt << std::endl;
    }
  }
}



