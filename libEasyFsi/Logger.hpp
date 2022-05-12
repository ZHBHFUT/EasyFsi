#pragma once

namespace EasyLib{

    void set_output_func(void (*func)(const char*));

    void error(const char* format, ...);
    void warn (const char* format, ...);
    void info (const char* format, ...);
    void debug(const char* format, ...);
}
