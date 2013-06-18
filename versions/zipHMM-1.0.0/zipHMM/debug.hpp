#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <string>
#include <cstdarg>

// #define _DEBUG 1

#ifdef _DEBUG 
  #define DEBUG(fmt, ...) debug(fmt, ##__VA_ARGS__)
#else
  #define DEBUG(fmt, ...)
#endif

void debug(const char * format, ...);

#endif
