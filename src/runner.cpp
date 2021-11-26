// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#include <cstdio>
#include <thread>
#include "runner.h"
#include "CountProfile.h"
#include "filehandling.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "LookuptableBase.h"
#include "util.h"
#include "Info.h"

#define SEQ_BUFSIZE 4096

KSEQ_INIT(int, read)

int processReads(string readsname,
                 LookupTableBase *lookuptable,
                 const KmerTranslator *translator,
                 int (*processCountProfile)(CountProfile &, void *),
                 void *processArgs,
                 int skip,
                 FILE *skipReadFile,
                 bool silent,
                 size_t chunkStart,
                 size_t chunkEnd)
{

  if (!silent)
    Info(Info::INFO) << "process reads from file " << readsname << "...\n";
  FILE *reads = openFileOrDie(readsname, "r");

  //TODO: check what happens if one thread have problems with a file, kill whole process then
  int fd = fileno(reads);
  kseq_t *seq = kseq_init(fd);

  CountProfile countprofile(translator, lookuptable);

  fseek(reads, chunkStart, SEEK_SET);
  /* iterate over every single fasta/fastq entry  */
  while (kseq_read(seq) >= 0) {

    const size_t len = seq->seq.l;
    const char *seqName = seq->name.s;

    /* fill profile */
    SequenceInfo *seqinfo = new SequenceInfo{seq->name.s, seq->comment.l!=0 ? string(seq->comment.s) : string(""), seq->seq.s,
                                             seq->qual.s!=NULL ? string(seq->qual.s) : string(""), seq->qual.s!=NULL ? '@':'>'};

    unsigned int kmerSpan = translator->getSpan();
    if (len < skip + kmerSpan) {
      if (skipReadFile)
        sequenceInfo2FileEntry(seqinfo, skipReadFile);
      else
        Info(Info::WARNING) << "WARNING: sequence " << seqName << " is too short, it'll be skipped\n";
      continue;
    }

    countprofile.fill(seqinfo, len);
    /* use function pointer for what to do with profile */
    processCountProfile(countprofile, processArgs);

    delete seqinfo;
    if (lseek(fd, 0, SEEK_CUR) > chunkEnd)
      break;
  }
  kseq_destroy(seq);
  fclose(reads);
  if (!silent)
    Info(Info::INFO) << "...completed\n";
  return EXIT_SUCCESS;
}

//int concatenateThreadResults(string resultFilename, string tresultFilenamePrefix, int threadNum, bool deleteTresults)
//{
//  //concatenate result file from threads temp results
//  FILE *resultFile = openFileOrDie(resultFilename, "w");
//  int fd_resultFile = fileno(resultFile);
//  for(size_t threadID = 0; threadID < threadNum; threadID++)
//  {
//    std::cerr << "write temporary results of thread " << threadID << " to resultfile " << std::endl << std::flush;
//    size_t size;
//    char buf[KS_BUFSIZE]; //TODO:change to using same buffer as kseq, or change kseq
//    FILE *source = openFileOrDie(tresultFilenamePrefix + std::to_string(threadID), "r");
//    int fd_source = fileno(source);
//    while ((size = read(fd_source, buf, KS_BUFSIZE)) > 0)
//    {
//      if(write(fd_resultFile, buf, size)==-1)
//      {
//        fprintf(stderr, "Error when writing to file %s: \n%s\n",
//                resultFilename.c_str(), strerror(errno));
//        fclose(source);
//        fclose(resultFile);
//        exit(EXIT_FAILURE);//TODO
//      }
//    }
//    fclose(source);
//    if (deleteTresults)
//      remove((tresultFilenamePrefix+ std::to_string(threadID)).c_str());
//  }
//  std::cerr << "Finish resultfile " << resultFilename << std::endl << std::flush;
//  fclose(resultFile);
//}

//
//int processreadsParallel(string readsname,
//                           string resultFilename,
//                           const LookupTableBase* lookuptable,
//                           const KmerTranslator* translator,
//                           int (*processCountProfile)(CountProfile &, FILE*),
//                           int threadNum)
//{
//
// /* string outdir = string("tmp/");
//  _mkdir(outdir)*/
//
//  FILE *reads = openFileOrDie(readsname, "r");
//  int fd = fileno(reads);
//  lseek(fd,0,SEEK_END);
//  uint64_t filesize = lseek(fd,0,SEEK_CUR);
//
//  /* chunksize as filesize divided by number of threads and then forced
//     to a multiple of KS_BUFSIZE (= kseq.h buffersize) */
//  uint64_t chunksize = ((uint64_t) ((filesize/threadNum)/KS_BUFSIZE))*KS_BUFSIZE;
//  if (chunksize == 0)
//    chunksize=std::min(filesize,(uint64_t) KS_BUFSIZE);
//
//  /* calculate concrete chunkStart to make sure chunks are approximatly equal
//   * in terms of bytes but don't break up fasta/fastq entries. Afterwards call
//   * function processReads for every thread on the specific chunk */
//  char buffer[SEQ_BUFSIZE];
//  uint64_t chunkStart = 0, chunkNextStart=0;
//  bool foundNextStart = true;
//  std::thread **threads = new std::thread*[threadNum];
//  size_t threadID=0;
//  string tresultFilename = resultFilename + string(".temp");
//  for(; threadID < threadNum && foundNextStart; threadID++)
//  {
//    {
//      lseek(fd, chunkStart+chunksize, SEEK_SET);
//      chunkNextStart = chunkStart+chunksize;
//      foundNextStart = false;
//
//      while(!foundNextStart)
//      {
//        size_t numElem = read(fd, buffer, SEQ_BUFSIZE);
//        //TODO: errno case numElem=-1
//        for(size_t jdx = 0; numElem !=0 && jdx < numElem-1; jdx++)
//        {
//          if (buffer[jdx] == '\n' &&
//              (buffer[jdx+1] == '>' || buffer[jdx+1] == '@'))
//          {
//            chunkNextStart += jdx+1;
//            foundNextStart = true;
//            break;
//          }
//        }
//        if (numElem < SEQ_BUFSIZE) //eof
//          break;
//
//        /* if the buffer does not contain any start of a fasta/fastq entry,
//         * increase start by the whole buffer size and read next 4096 elements*/
//        if (!foundNextStart)
//          chunkNextStart += SEQ_BUFSIZE;
//      }
//    }
//
//    uint64_t chunkEnd = chunkStart+chunksize;
//
//    /* special case: chunkEnd of last thread have to be the fileEnd */
//    if(threadID == threadNum-1)
//      chunkEnd = filesize;
//
//
//    std::cerr << "thread " << threadID << " working on [" << chunkStart << ","
//              << chunkStart+chunksize << "]" << std::endl << std::flush;
//
//    threads[threadID] = new std::thread(processReads,
//                                   readsname,
//                                   tresultFilename+std::to_string(threadID),//TODO:store name in dic?,
//                                   lookuptable,
//                                   translator,
//                                   processCountProfile,
//                                   chunkStart,
//                                   chunkEnd);
//
//    chunkStart = chunkNextStart;
//  }
//
//  /* join and delete threads */
//  threadNum=threadID;
//  for(threadID = 0; threadID < threadNum; threadID++)
//  {
//    threads[threadID]->join();
//    delete threads[threadID];
//    std::cerr << "thread " << threadID << " finished" << std::endl << std::flush;
//  }
//  delete[] threads;
//  fclose(reads);
//
//  /* collate thread iterim results */
//  concatenateThreadResults(resultFilename, tresultFilename, threadNum, true);
//}
