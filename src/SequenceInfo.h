//Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#ifndef SEQUENCEINFO_H
#define SEQUENCEINFO_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

enum outMode {AUTO, FASTA};

using namespace std;
typedef struct {
  string name, comment, seq, qual;
  char sep;
}
SequenceInfo;



inline void sequenceInfo2FileEntry(SequenceInfo *seqinfo, FILE *fp, outMode mode=AUTO) {


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
