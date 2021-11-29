//Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#ifndef SEQUENCEINFO_H
#define SEQUENCEINFO_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "Info.h"
#include "util.h"
#include "filehandling.h"

enum SeqInfoMode {AUTO, FASTA, FASTQ};

using namespace std;
typedef struct {
  string name, comment, seq, qual;
  char sep;
}
SequenceInfo;

inline SeqInfoMode getSeqMode(std::string filename) {
  FILE *fp = openFileOrDie(filename, "r");
  char c;
  if (fscanf(fp,"%c", &c) != 1) {
    Info(Info::ERROR) << "ERROR: Can not read from " << filename << "\n";
    EXIT(EXIT_FAILURE);
  }
  SeqInfoMode mode;
  if (c == '>')
    mode = FASTA;
  else if (c == '@')
    mode = FASTQ;
  else {
    Info(Info::ERROR)  << "ERROR: Invalid first character found in " << filename << "expect '>' or '@'\n";
    EXIT(EXIT_FAILURE);
  }

  fclose(fp);
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

#endif //SEQUENCEINFO_H
