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

    string tresultFilename = resultFilename + string(".temp") + std::to_string(threadID);//TODO:store name in dic?
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
  
  threadNum=threadID;
  for(threadID = 0; threadID < threadNum; threadID++)
  {
    threads[threadID]->join();
    delete threads[threadID];
    std::cerr << "thread " << threadID << "finished" << std::endl << std::flush;
  }
  delete[] threads;
  fclose(seqFile);

  //TODO: fopen or open, fread/fwrite or read/write
  //concatenate result file from threads temp results
  FILE *resultFile = openFileOrDie(resultFilename, "w");
  for(threadID = 0; threadID < threadNum; threadID++)
  {
    std::cerr << "write temporary results of thread " << threadID << " to resultfile " << std::endl << std::flush;
    size_t size;
    char buf[BUFSIZ]; //TODO:change to using same buffer as above, or change above
    FILE *source = openFileOrDie(resultFilename + string(".temp") + std::to_string(threadID), "r");
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



int processSeqFile(string seqFilename,
                   string resultFilename,
                   const Lookuptable *lookuptable,
                   const KmerTranslator *translator,
                   int (*processCountProfile)(CountProfile &, FILE*),
                   size_t chunkStart,
                   size_t chunkEnd)
{

  FILE *seqFile = openFileOrDie(seqFilename, "r");//TODO: lower level?
  //TODO: check what happens if one thread have problems wit a file, kill whole process then
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



