#ifndef TOOL
#define TOOL

struct ToolInfo {
    const char *cmd;
    int (*callerFunction)(int, const char **, const struct ToolInfo*);
    //std::vector<cocoParameter>* params;
    const char *decript;
    const char *author;
    const char *usage;
};

#endif // TOOL

