#include <thread>
#include "processSequences.h"
#include "CountProfile.h"
#include "filehandling.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "Lookuptable.h"
#include "util.h"

KSEQ_INIT(int, read)

int testpara(size_t id, const char* seqFilename, size_t chunkStart, size_t chunkEnd)
{
  FILE *seqFile = openFileOrDie(seqFilename, "r");
  fseek(seqFile, chunkStart, SEEK_SET);
  int fd = fileno(seqFile);
  kseq_t *seq = kseq_init(fileno(seqFile));
  FILE *resultFile = openFileOrDie(std::to_string(id), "w");

  /* iterate over every single read  */
  while (kseq_read(seq) >= 0)
  {

    fprintf(resultFile, "@%s %s\n",seq->name.s,seq->comment.s);
    fprintf(resultFile, "%s\n",seq->seq.s);
    fprintf(resultFile, "+\n%s\n", seq->qual.s);
    //fprintf(resultFile, "%u\n", lseek(fd,0,SEEK_CUR));
    if (lseek(fd,0,SEEK_CUR) > chunkEnd)
      break;
  }

  kseq_destroy(seq);
  fclose(seqFile);
  fclose(resultFile);
}

int processSeqFileParallel(string seqFilename,
                           string resultFilename,
                           const Lookuptable* lookuptable,
                           const KmerTranslator* translator,
                           int (*processCountProfile)(CountProfile &, FILE*),
                           int threadNum)
{

 /* string outdir = string("tmp/");
  _mkdir(outdir)*/

  FILE *seqFile = openFileOrDie(seqFilename, "r");
  fseek(seqFile,0,SEEK_END);
  uint64_t filesize = ftell(seqFile);
  char buffer[4096]; //TODO: macro, increase

  // chunksize as filesize divided by number of threads and then forced
  // to a multiple of 4096 (= kseq buffersize)
  uint64_t chunksize = (filesize/threadNum >> 12) << 12;
  if (chunksize == 0)
    chunksize=std::min(filesize,(uint64_t) 4096);

  //TODO: special case chunksize 0

  std::thread **threads = new std::thread*[threadNum];
  uint64_t chunkStart = 0, chunkNextStart=0;
  bool foundNextStart = true;
  size_t threadID=0;
  for(; threadID < threadNum && foundNextStart; threadID++)
  {
    //chunkEnd = chunkStart+chunksize;
    //if (chunkEnd + 4096 < filesize)
    {
      fseek(seqFile, chunkStart+chunksize, SEEK_SET);
      chunkNextStart = chunkStart+chunksize;
      foundNextStart = false;

      while(!foundNextStart)
      {
        size_t numElem = fread(buffer, sizeof(char), 4096, seqFile);
        for(size_t jdx = 0; numElem !=0 && jdx < numElem-1; jdx++)
        {
          if (buffer[jdx] == '\n' &&
              (buffer[jdx+1] == '>' || buffer[jdx+1] == '@'))
          {
            chunkNextStart += jdx+1;
            foundNextStart = true;
            break;
          }
        }
        if (feof(seqFile))
          break;

        if (!foundNextStart)
          chunkNextStart += 4096;
      }
    }

    uint64_t chunkEnd = chunkStart+chunksize;
    if(threadID == threadNum-1)
      chunkEnd=filesize;


    std::cerr << "thread " << threadID << " working on [" << chunkStart << ","
              << chunkStart+chunksize << "]" << std::endl << std::flush;

    string tresultFilename = resultFilename + string(".temp") + std::to_string(threadID);
    threads[threadID] = new std::thread(processSeqFile,
                                   seqFilename,
                                   tresultFilename,
                                   lookuptable,
                                   translator,
                                   processCountProfile,
                                   chunkStart,
                                   chunkEnd);

    chunkStart = chunkNextStart;
  }
  for(int i = 0; i < threadID; i++)
  {
    threads[i]->join();
    delete threads[i];
  }
  delete[] threads;

  //concatenate result files from threads temp results
}



int processSeqFile(string seqFilename,
                   string resultFilename,
                   const Lookuptable *lookuptable,
                   const KmerTranslator *translator,
                   int (*processCountProfile)(CountProfile &, FILE*),
                   size_t chunkStart,
                   size_t chunkEnd)
{

  FILE *seqFile = openFileOrDie(seqFilename, "r");
  int fd = fileno(seqFile);
  kseq_t *seq = kseq_init(fd);
  FILE *resultFile = openFileOrDie(resultFilename, "w");
  CountProfile countprofile(translator, lookuptable);

  fseek(seqFile, chunkStart, SEEK_SET);

  /* iterate over every single read  */
  while (kseq_read(seq) >= 0)
  {
    const size_t len = seq->seq.l;
    const char* readSeq = seq->seq.s;
    const char* readName = seq->name.s;
    SeqType seqStr; seqStr.reserve(len);

    unsigned int kmerSpan = translator->getSpan();
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

    /* use function pointer for what to do with countprofile */
    processCountProfile(countprofile, resultFile);

    if (lseek(fd,0,SEEK_CUR) > chunkEnd)
       break;

  }
  kseq_destroy(seq);
  fclose(seqFile);
  fclose(resultFile);
  return EXIT_SUCCESS;
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
                                 lookuptable,
                                 translator,
                                 processCountProfile);
    if (retval != EXIT_SUCCESS)
    {
      std::cerr << "ERROR processing read file " << *sampleIt << std::endl;
    }
  }
}



