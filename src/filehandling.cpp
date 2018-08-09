#include "filehandling.h"

void _mkdir(const char *dir, mode_t mode)
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
void _mkdir(string &dir, mode_t)
{
  std::cout << "creating dir " << dir << std::endl << std::flush;
  _mkdir(dir.c_str());
}

vector<string> *getFileList(const char *fileListFilename)
{
  vector<string> *fileList = new vector<string>();
  string line;
  ifstream fp;
  fp.open(fileListFilename);
  while (getline (fp,line))
    fileList->push_back(line);
  fp.close();

  return fileList;
}


