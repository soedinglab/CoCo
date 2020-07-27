//Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>

#ifndef UTIL
#define UTIL

//#define DEBUG

#ifndef EXIT
#define EXIT(exitCode) do { int __status = (exitCode);\
                            std::cerr.flush(); std::cout.flush();\
                            exit(__status); } while(0)
#endif

enum UsageMode {
  SIMPLE = 0,
  EXTENDED
};

#endif // UTIL

