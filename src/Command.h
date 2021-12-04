// Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>
#ifndef COMMAND
#define COMMAND

#include <vector>

enum coco_commands {
  CORRECTOR = 1, CONSENSUS = 2, FILTER = 4, ABUNDANCE_ESTIMATOR = 8, PROFILE = 16, COUNTS2FLAT = 32
};
struct cocoOption;

typedef struct Command {
  const char *cmd;

  int (*callerFunction)(int, const char **, const struct Command *);

  std::vector<cocoOption*> *opt;
  const char *descriptShort;
  const char *descriptLong;
  const char *author;
  const char *usage;
  const char id;
} Command;

#endif // COMMAND

