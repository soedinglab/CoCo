#ifndef UTIL
#define UTIL

#ifndef EXIT
#define EXIT(exitCode) do { int __status = (exitCode);\
                            std::cerr.flush(); std::cout.flush();\
                            exit(__status); } while(0)
#endif

#endif // UTIL

