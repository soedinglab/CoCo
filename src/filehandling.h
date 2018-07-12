#ifndef FILEHANDLING_H
#define FILEHANDLING_H
#include <sys/stat.h>
#include <string>



void _mkdir(const char *dir, mode_t mode=S_IRWXU | S_IRWXG | S_IROTH)
{
  char tmp[256];
  char *p = NULL;
  size_t len;

  snprintf(tmp, sizeof(tmp),"%s",dir);
  len = strlen(tmp);
  if(tmp[len - 1] == '/')
    tmp[len - 1] = 0;
  for(p = tmp + 1; *p; p++)
  {
    if(*p == '/')
    {
      *p = 0;
      mkdir(tmp, mode);
      *p = '/';
    }
  }
  mkdir(tmp, mode);
}
void _mkdir(string &dir, mode_t= S_IRWXU | S_IRWXG | S_IROTH)
{
  std::cout << "creating dir " << dir << std::endl << std::flush;
  _mkdir(dir.c_str());
}
#endif // FILEHANDLING_H
