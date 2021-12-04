// Written by Annika Jochheim <annika.jochheim@mpibpc.mpg.de>
#ifndef FILEHANDLING_H
#define FILEHANDLING_H

#include <sys/stat.h>
#include <string>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

bool _mkdir(const char *dir, mode_t mode = S_IRWXU | S_IRWXG | S_IROTH);

bool _mkdir(std::string &dir, mode_t= S_IRWXU | S_IRWXG | S_IROTH);

bool fileExists(const char* fileName);

bool directoryExists(const char* dirName);

FILE *openFileOrDie(const std::string &fileName, const char *mode);

FILE *openFileOrDie(const char *fileName, const char *mode);


string getFilename(string wholeFilePath);

string getFileExtension(string wholeFilePath);

vector<std::string> *getFileList(const char *fileListFilename);


#endif // FILEHANDLING_H
