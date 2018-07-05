#ifndef TOOL
#define TOOL

#include <vector>

struct cocoOption;

typedef struct ToolInfo {
    const char *cmd;
    int (*callerFunction)(int, const char **, const struct ToolInfo*);
    std::vector<cocoOption>* opt;
    const char *descriptShort;
    const char *descriptLong;
    const char *author;
    const char *usage;
    unsigned int toolNum;
} ToolInfo;

#endif // TOOL

