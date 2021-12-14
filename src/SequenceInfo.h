//Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>

#ifndef SEQUENCEINFO_H
#define SEQUENCEINFO_H

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "Info.h"
#include "util.h"
#include "filehandling.h"
#include "KSeqWrapper.h"
enum SeqInfoMode {AUTO, FASTA, FASTQ};

using namespace std;
typedef struct {
  string name, comment, seq, qual;
  char sep;
}
SequenceInfo;

inline SeqInfoMode getSeqMode(std::string filename) {

  KSeqWrapper *kseqReader = KSeqFactory(filename.c_str());


  if (!kseqReader->ReadEntry()) {
    Info(Info::ERROR) << "ERROR: Can not read from " << filename << "\n";
    EXIT(EXIT_FAILURE);
  }
  const KSeqWrapper::KSeqEntry &seq = kseqReader->entry;
  SeqInfoMode mode;
  if (seq.qual.l == 0)
    mode = FASTA;
  else if (seq.qual.l == seq.sequence.l)
    mode = FASTQ;
  else {
    Info(Info::ERROR)  << "ERROR: Invalid first entry found in " << filename << "\n";
    EXIT(EXIT_FAILURE);
  }

  delete kseqReader;
  return mode;
}

inline void sequenceInfo2FileEntry(SequenceInfo *seqinfo, FILE *fp, SeqInfoMode mode=AUTO) {


  fwrite(mode == AUTO ? &seqinfo->sep : ">", sizeof(char), 1, fp);
  fwrite(seqinfo->name.c_str(), sizeof(char), seqinfo->name.size(), fp);
  if (seqinfo->comment.size() > 0) {
    fwrite(" ", sizeof(char), 1, fp);
    fwrite(seqinfo->comment.c_str(), sizeof(char), seqinfo->comment.size(), fp);
  }
  fwrite("\n", sizeof(char), 1, fp);
  fwrite(seqinfo->seq.c_str(), sizeof(char), seqinfo->seq.size(), fp);
  fwrite("\n", sizeof(char), 1, fp);
  if (mode == AUTO && seqinfo->sep == '@') {
    fwrite("+\n", sizeof(char), 2, fp);
    fwrite(seqinfo->qual.c_str(), sizeof(char), seqinfo->qual.size(), fp);
    fwrite("\n", sizeof(char), 1, fp);
  }
}

inline char getAvgQual(std::string qual, unsigned int pos) {

  int qval = 0, count = 0;
  if(!qual.empty()) {
    if (pos > 0) {
      qval += qual[pos - 1];
      count++;
    }
    if (pos < qual.length()-1) {
      qval += qual[pos + 1];
      count++;
    }
    return (char) (qval/count);
  }

  return 33; //error value, phred 33
}

#endif //SEQUENCEINFO_H
