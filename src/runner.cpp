// Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>

#include <cstdio>
#include <thread>
#include "runner.h"
#include "CountProfile.h"
#include "filehandling.h"
#include "KmerTranslator.h"
#include "kseq.h"
#include "LookuptableBase.h"
#include "Info.h"
#include "KSeqWrapper.h"
//#define SEQ_BUFSIZE 4096


int processReads(const string &readsname,
                 LookupTableBase *lookuptable,
                 const KmerTranslator *translator,
                 int (*processCountProfile)(CountProfile &, void *, bool),
                 void *processArgs,
                 int skip,
                 bool silent)
{

  if (!silent)
    Info(Info::INFO) << "...process reads from file " << readsname << "\n";


  KSeqWrapper* kseqReader = KSeqFactory(readsname.c_str());
  CountProfile countprofile(translator, lookuptable);

  /* iterate over every single fasta/fastq entry  */
  while(kseqReader->ReadEntry()) {
    const KSeqWrapper::KSeqEntry &seq = kseqReader->entry;
    if (seq.sequence.l == 0) {
      Info(Info::ERROR) << "ERROR: Invalid sequence record found\n";
      return EXIT_FAILURE;
    }

    /* fill profile */
    SequenceInfo *seqinfo = new SequenceInfo{seq.name.s,
                                    seq.comment.l!=0 ? string(seq.comment.s) : string(""),
                                             seq.sequence.s,seq.qual.s!=NULL ? string(seq.qual.s) : string(""),
                                        seq.qual.s!=NULL ? '@':'>'};

    unsigned int kmerSpan = translator->getSpan();
    if (seq.sequence.l < skip + kmerSpan) {

      countprofile.setSeqInfo(seqinfo);
      /* use function pointer for to skip sequence */
      processCountProfile(countprofile, processArgs, true);
    } else {

      countprofile.fill(seqinfo);
      /* use function pointer for what to do with profile */
      processCountProfile(countprofile, processArgs, false);
    }

    delete seqinfo;
  }
  delete kseqReader;
  if (!silent)
    Info(Info::INFO) << "...completed\n";
  return EXIT_SUCCESS;
}

int processPairedReads(const string &forwardReads,
                       const string &reverseReads,
                       LookupTableBase *lookuptable,
                       const KmerTranslator *translator,
                       int (*processCountProfile)(CountProfile &, CountProfile &,void *, bool),
                       void *processArgs,
                       int skip)
{

  Info(Info::INFO) << "...process paired reads from files " << forwardReads << " and " << reverseReads << "\n";

  KSeqWrapper* r1_kseq = KSeqFactory(forwardReads.c_str());
  KSeqWrapper* r2_kseq = KSeqFactory(reverseReads.c_str());

  CountProfile r1_countprofile(translator, lookuptable);
  CountProfile r2_countprofile(translator, lookuptable);

  int res1=0, res2=0, count=0;
  KSeqWrapper::KSeqEntry *r1;
  KSeqWrapper::KSeqEntry *r2;
  /* iterate over every single fasta/fastq entry  */
  while (r1_kseq->ReadEntry() && r2_kseq->ReadEntry()) {

    const KSeqWrapper::KSeqEntry &r1 = r1_kseq->entry;
    const KSeqWrapper::KSeqEntry &r2 = r2_kseq->entry;
    count++;

    const size_t len1 = r1.sequence.l;

    if (r1.sequence.l == 0 || r2.sequence.l == 0) {
      Info(Info::ERROR) << "ERROR: Invalid sequence record found when processing paired reads\n";
      return EXIT_FAILURE;
    }

    /* fill profile */
    SequenceInfo *r1_seqinfo = new SequenceInfo{r1.name.s, r1.comment.l!=0 ? string(r1.comment.s) : string(""), r1.sequence.s,
                                                r1.qual.s!=NULL ? string(r1.qual.s) : string(""), r1.qual.s!=NULL ? '@':'>'};
    SequenceInfo *r2_seqinfo = new SequenceInfo{r2.name.s, r2.comment.l!=0 ? string(r2.comment.s) : string(""), r2.sequence.s,
                                                r2.qual.s!=NULL ? string(r2.qual.s) : string(""), r2.qual.s!=NULL ? '@':'>'};


    unsigned int kmerSpan = translator->getSpan();
    if (r1.sequence.l < skip + kmerSpan || r2.sequence.l < skip + kmerSpan) {

      r1_countprofile.setSeqInfo(r1_seqinfo);
      r2_countprofile.setSeqInfo(r2_seqinfo);

      /* use function pointer for to skip sequence */
      processCountProfile(r1_countprofile, r2_countprofile, processArgs, true);
    } else {

      r1_countprofile.fill(r1_seqinfo);
      r2_countprofile.fill(r2_seqinfo);
      /* use function pointer for what to do with profile */
      processCountProfile(r1_countprofile, r2_countprofile, processArgs, false);
    }

    delete r1_seqinfo;
    delete r2_seqinfo;

  }

  delete r1_kseq;
  delete r2_kseq;

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
