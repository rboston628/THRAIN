#ifndef STRINGFORMATH
#define STRINGFORMATH

#include <string>
#include <cstdio>
#include <cstdarg>

std::string addstring (const std::string& s, const char* c);

std::string addstring (const char* c, const std::string& s);

std::string addstring (const char* c, const char* x);

std::string strmakef(const char *fmt, ...);

#endif