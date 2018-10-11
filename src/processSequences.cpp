#include <thread>
#include "processSequences.h"
#include "CountProfile.h"
#include "filehandling.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "Lookuptable.h"
#include "util.h"

KSEQ_INIT(int, read)

int processSeqFile(string seqFilename,
                   string resultFilename,
                   const Lookuptable *lookuptable,
                   const KmerTranslator *translator,
                   int (*processCountProfile)(CountProfile &, FILE*),
                   size_t chunkStart,
                   size_t chunkEnd) //TODO: check default
{

  FILE *seqFile = openFileOrDie(seqFilename, "r");//TODO: lower level?
  //TODO: check what happens if one thread have problems wit a file, kill whole process then
  int fd = fileno(seqFile);
  kseq_t *seq = kseq_init(fd);
  FILE *resultFile = openFileOrDie(resultFilename, "w");
  CountProfile countprofile(translator, lookuptable);

  fseek(seqFile, chunkStart, SEEK_SET);

  /* iterate over every single fasta/fastq entry  */
  while (kseq_read(seq) >= 0)
  {
    const size_t len = seq->seq.l;
    const char* seqNuc = seq->seq.s;
    const char* seqName = seq->name.s;
    SeqType seqStr; seqStr.reserve(len);

    unsigned int kmerSpan = translator->getSpan();
    if(len < kmerSpan)
    {
      fprintf(stderr, "WARNING: sequence %s is too short, it'll be skipped\n", seqName);
      continue;
    }

    /* sequence to 2bit representation */
    for(unsigned int pos = 0; pos < len; pos++)
    {
      seqStr.push_back((char)res2int[(int)seqNuc[pos]]);
    }

    /* fill profile */
    countprofile.fill(seqStr, seqName);

    /* use function pointer for what to do with profile */
    processCountProfile(countprofile, resultFile);

    if (lseek(fd,0,SEEK_CUR) > chunkEnd)
       break;

  }
  kseq_destroy(seq);
  fclose(seqFile);
  fclose(resultFile);
  return EXIT_SUCCESS;
}

int concatenateThreadResults(string resultFilename, string tresultFilenamePrefix, int threadNum)
{
  //TODO: fopen or open, fread/fwrite or read/write
  //concatenate result file from threads temp results
  FILE *resultFile = openFileOrDie(resultFilename, "w");
  for(size_t threadID = 0; threadID < threadNum; threadID++)
  {
    std::cerr << "write temporary results of thread " << threadID << " to resultfile " << std::endl << std::flush;
    size_t size;
    char buf[BUFSIZ]; //TODO:change to using same buffer as kseq, or change kseq
    FILE *source = openFileOrDie(tresultFilenamePrefix + std::to_string(threadID), "r");
    while ((size = fread(buf, sizeof(char), BUFSIZ, source)) > 0)
    {
      fwrite(buf, sizeof(char), size, resultFile);
      if(ferror(resultFile))
      {
        fprintf(stderr, "Error when writing to file %s: \n%s\n",
                resultFilename.c_str(), strerror(errno));
        fclose(source);
        fclose(resultFile);
        exit(EXIT_FAILURE);//TODO
      }
    }
    fclose(source);
  }
  std::cerr << "Finish resultfile " << resultFilename << std::endl << std::flush;
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

  /* chunksize as filesize divided by number of threads and then forced
     to a multiple of 4096 (= kseq buffersize) */
  uint64_t chunksize = (filesize/threadNum >> 12) << 12;
  if (chunksize == 0)
    chunksize=std::min(filesize,(uint64_t) 4096);

  /* calculate concrete chunkStart to make sure chunks are approximatly equal
   * in terms of bytes but don't break up fasta/fastq entries. Afterwards call
   * function processSeqFile for every thread on the specific chunk */
  std::thread **threads = new std::thread*[threadNum];
  uint64_t chunkStart = 0, chunkNextStart=0;
  bool foundNextStart = true;
  size_t threadID=0;
  string tresultFilename = resultFilename + string(".temp");
  for(; threadID < threadNum && foundNextStart; threadID++)
  {
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

        /* if the buffer does not contain any start of a fasta/fastq entry,
         * increase start by the whole buffer size and read next 4096 elements*/
        if (!foundNextStart)
          chunkNextStart += 4096;
      }
    }

    uint64_t chunkEnd = chunkStart+chunksize;

    /* special case: chunkEnd of last thread have to be the fileEnd */
    if(threadID == threadNum-1)
      chunkEnd = filesize;


    std::cerr << "thread " << threadID << " working on [" << chunkStart << ","
              << chunkStart+chunksize << "]" << std::endl << std::flush;

    threads[threadID] = new std::thread(processSeqFile,
                                   seqFilename,
                                   tresultFilename+std::to_string(threadID),//TODO:store name in dic?,
                                   lookuptable,
                                   translator,
                                   processCountProfile,
                                   chunkStart,
                                   chunkEnd);

    chunkStart = chunkNextStart;
  }

  /* join and delete threads */
  threadNum=threadID;
  for(threadID = 0; threadID < threadNum; threadID++)
  {
    threads[threadID]->join();
    delete threads[threadID];
    std::cerr << "thread " << threadID << " finished" << std::endl << std::flush;
  }
  delete[] threads;
  fclose(seqFile);

  /* collate thread iterim results */
  concatenateThreadResults(resultFilename, tresultFilename, threadNum);
}
