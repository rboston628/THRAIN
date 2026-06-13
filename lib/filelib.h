// methods for interacting with the file system

#ifndef THRAINFILELIB
#define THRAINFILELIB

#include <string>

namespace filelib {
void makedir(std::string const &path);
void remove(std::string const &path);
void touch(std::string const &path);
bool exists(std::string const &path);
}

#endif
