// inspired by solutions on 
// https://codereview.stackexchange.com/questions/187183/create-a-c-string-using-printf-style-formatting

#ifndef STRINGFORMATCPP
#define STRINGFORMATCPP

#include "string.h"
#include <cstdio>
#include <cstdarg>

std::string addstring (const std::string& s, const char* c){
    std::string r = s;
    for(std::size_t i=0; c[i]!=0; i++){
        r += c[i];
    }
    return r;
}

std::string addstring (const char* c, const std::string& s){
    std::string r = "";
    for(std::size_t i=0; c[i]!=0; i++){
        r += c[i];
    }
    r += s;
    return r;
}

std::string addstring (const char* c, const char* x){
    std::string r = "";
    for(std::size_t i=0; c[i]!=0; i++){
        r += c[i];
    }
    for(std::size_t j=0; x[j]!=0; j++){
        r += x[j];
    }
    return r;
}

std::string strmakef(const char *fmt, ...)
{
    char buf[256];

    va_list args;
    va_start(args, fmt);
    const auto r = std::vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);

    if (r < 0)
        // conversion failed
        return {};

    const size_t len = r;
    if (len < sizeof buf)
        // we fit in the buffer
        return { buf, len };

    std::string s(len, '\0');
    va_start(args, fmt);
    std::vsnprintf(&(*s.begin()), len+1, fmt, args);
    va_end(args);
    return s;
}



#endif