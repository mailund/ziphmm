#include "debug.hpp"

#include <cstdarg>
#include <stdio.h>

void debug(const char * format, ...) {
    char buffer[256];
    va_list args;
    va_start (args, format);
    vsprintf (buffer, format, args);
    printf("%s", buffer);
    va_end (args);
}
