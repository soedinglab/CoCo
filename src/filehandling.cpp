// Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>
#include "filehandling.h"
#include "util.h"
#include "Info.h"

bool _mkdir(const char *dir, mode_t mode) {
  char tmp[1024];
  char *p = NULL;
  size_t len;

  snprintf(tmp, sizeof(tmp), "%s", dir);
  len = strlen(tmp);
  if (tmp[len - 1] == '/')
    tmp[len - 1] = 0;
  for (p = tmp + 1; *p; p++) {
    if (*p == '/') {
      *p = 0;
      if(mkdir(tmp, mode) != 0) {
        return false;
      }
      *p = '/';
    }
  }
  if(mkdir(tmp, mode) != 0) {
    return false;
  }
  return true;
}

bool _mkdir(string &dir, mode_t) {
  Info(Info::INFO) << "Creating dir " << dir << "\n";
  return _mkdir(dir.c_str());
}

bool fileExists(const char* fileName) {
  struct stat st;
  return stat(fileName, &st) == 0;
}

bool directoryExists(const char* dirName) {
  struct stat st;
  return stat(dirName, &st) == 0 && S_ISDIR(st.st_mode);
}

FILE *openFileOrDie(std::string fileName, const char *mode) {
  return openFileOrDie(fileName.c_str(), mode);
}

FILE *openFileOrDie(const char *fileName, const char *mode) {

  FILE *file;
  file = fopen(fileName, mode);
  if (file == NULL) {
    Info(Info::ERROR) << "ERROR: opening failed for file " << fileName << "\n";
    perror(fileName);
    EXIT(EXIT_FAILURE);
  }
  return file;
}

string getFileExtension(string wholeFilePath) {
  string ext = "";
  size_t lastdot = wholeFilePath.find_last_of(".");
  if (lastdot != std::string::npos)
    ext = wholeFilePath.substr(lastdot, wholeFilePath.size());
  return (ext);
}

string getFilename(string wholeFilePath) {
  string filename = wholeFilePath;
  size_t lastdot = filename.find_last_of(".");
  if (lastdot != std::string::npos)
    filename = filename.substr(0, lastdot);

  return string(basename(filename.c_str()));
}

vector<string> *getFileList(const char *fileListFilename) {
  vector<string> *fileList = new vector<string>();
  string line;
  ifstream fp;
  fp.open(fileListFilename);
  while (getline(fp, line))
    fileList->push_back(line);
  fp.close();

  return fileList;
}


