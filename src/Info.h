//
// Written by Annika Seidel <annika.seidel@mpibpc.mpg.de>
//

#include <string>
#include <iostream>

#ifndef INFO_H
#define INFO_H

class Info {

private:
  const int level;
  std::string buffer;

public:
  static const int ERROR = 1;
  static const int WARNING = 2;
  static const int INFO = 3;
  static const int DEBUG = 4;

  static int verboseLevel;

  explicit Info( int level ) : level(level) {};

  ~Info() {};

  static void setVerboseLevel (int i);

  template<typename T>
  Info& operator<<(const T& t)
  {
    if (level <= verboseLevel) {
      if (level == INFO) {
        std::cout << t;
      }
      else
        std::cerr << t;
    }
    return *this;
  }
};



#endif //INFO_H
