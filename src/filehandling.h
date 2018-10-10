#ifndef FILEHANDLING_H
#define FILEHANDLING_H
#include <sys/stat.h>
#include <string>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

void _mkdir(const char *dir, mode_t mode=S_IRWXU | S_IRWXG | S_IROTH);

void _mkdir(std::string &dir, mode_t= S_IRWXU | S_IRWXG | S_IROTH);

FILE* openFileOrDie(std::string fileName, const char * mode);

FILE* openFileOrDie(const char *fileName, const char * mode);


string get_filename(string wholeFilePath);

vector<std::string> *getFileList(const char *fileListFilename);


#endif // FILEHANDLING_H
