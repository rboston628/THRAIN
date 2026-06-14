#include "filelib.h"

#include <cerrno>
#include <sys/stat.h>
#if defined(_WIN32)
#include <direct.h>
#endif

#include <stdexcept>
#include <cstdio>
#include <string>
#include <algorithm>

namespace {
std::string _normal_to_windows(std::string const &path) {
  std::string windows_path = path;
  std::replace(windows_path.begin(), windows_path.end(), '/', '\\');
  return windows_path;
}
}

namespace filelib {
void makedir(std::string const &path) {
  // get the parent, recurse baxkward until parent path exists, then create forward
  std::string const parent = path.substr(0, path.find_last_of("/\\"));
  if (!parent.empty() && !exists(parent)) {
    makedir(parent);
  }
  // make the final directory
#if defined(_WIN32)
  std::string const windows_path = _normal_to_windows(path);
  if (_mkdir(windows_path.c_str()) != 0 && errno != EEXIST) throw std::runtime_error("mkdir failed to create directory " + std::string(path));
#else
  if (mkdir(path.c_str(), 0777) != 0 && errno != EEXIST) throw std::runtime_error("mkdir failed to create directory " + std::string(path));
#endif
}

void remove(std::string const &path) {
  int result;
  if (!exists(path)) return;
#if defined(_WIN32)
  std::string const windows_path = _normal_to_windows(path);
  if (0 != (result = std::remove(windows_path.c_str()))) {
    result = system(("rmdir /s /q " + windows_path).c_str());
  }
#else
  if (0 != (result = std::remove(path.c_str()))) {
    result = system(("rm -rf " + path).c_str());
  }
#endif
  if (result != 0) throw std::runtime_error("Failed to remove " + path);
}

void touch(std::string const &path) {
  FILE *fp = fopen(path.c_str(), "a");
  if (!fp) throw std::runtime_error("Failed to touch file " + path);
  fclose(fp);
}

bool exists(std::string const &path) {
  struct stat buffer;
  return (stat(path.c_str(), &buffer) == 0);
}
}