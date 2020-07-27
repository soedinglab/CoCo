// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
#ifndef COMMAND
#define COMMAND


#include <vector>

enum coco_commands {
  PROFILE = 1, FILTER = 2, ABUNDANCE_ESTIMATOR = 8, CONSENSUS = 16
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

