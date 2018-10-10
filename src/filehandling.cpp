#include "filehandling.h"
#include "util.h"

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


FILE* openFileOrDie(std::string fileName, const char * mode)
{
  return openFileOrDie(fileName.c_str(), mode);
}

FILE* openFileOrDie(const char *fileName, const char * mode)
{

  FILE* file;
  file = fopen(fileName, mode);
  if(file == NULL)
  {
    fprintf(stderr, "ERROR: opening failed for file %s\n", fileName);
    perror(fileName);
    EXIT(EXIT_FAILURE);
  }
  return file;
}

string get_filename(string wholeFilePath)
{
  string filename = wholeFilePath;
  size_t lastdot = filename.find_last_of(".");
  if (lastdot != std::string::npos)
    filename = filename.substr(0, lastdot);

  return string(basename(filename.c_str()));
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


