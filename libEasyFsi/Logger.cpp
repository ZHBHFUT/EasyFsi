/*-------------------------------------------------------------------------
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from
the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not
   be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution.
-------------------------------------------------------------------------*/

//!-------------------------------------------------------------
//! @file       Logger.cpp
//!             The implement of logging functions.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------
#include "Logger.hpp"

#if defined _WIN32
#pragma message ("platform: windows")
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>  // see enable_console_color
#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#endif
#endif

#include <cstdio>
#include <cstdarg>
#include <string>
#include <stdexcept>

#include "Assert.hpp"

namespace EasyLib {
    static void (*print)(const char*) = nullptr;
    static void (*fexit)() = nullptr;

#if defined _WIN32
    static void enable_console_color_for_windows()
    {
        // References:
        //SetConsoleMode() and ENABLE_VIRTUAL_TERMINAL_PROCESSING?
        //https://stackoverflow.com/questions/38772468/setconsolemode-and-enable-virtual-terminal-processing

        // Windows console with ANSI colors handling
        // https://superuser.com/questions/413073/windows-console-with-ansi-colors-handling

        static bool init = false;
        if (init)return;

        HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
        DWORD dwMode = 0;
        GetConsoleMode(hOut, &dwMode);
        dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
        SetConsoleMode(hOut, dwMode);

        HANDLE hErr = GetStdHandle(STD_ERROR_HANDLE);
        dwMode = 0;
        GetConsoleMode(hErr, &dwMode);
        dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
        SetConsoleMode(hErr, dwMode);

        init = true;
    }
#endif
    static std::string format_(const char* format, va_list args, const char* prefix = nullptr, const char* postfix = nullptr)
    {
        std::string output((prefix && *prefix) ? prefix : "");
        
#ifdef _MSC_VER
        //use _vscprintf, vsprintf_s
        auto msg_len = _vscprintf(format, args);// not counting the terminating null character.
        if (msg_len <= 0)return output;

        auto offset = output.length();
        output.resize(offset + (size_t)msg_len);// string[string.size()] is guaranteed to be '\0'
        //vsprintf(output.data(), fmt, arg_ptr);//C4996
        vsprintf_s(&output[offset], msg_len + 1, format, args);
#else
        //use vsnprintf, vsprintf
        auto msg_len = vsnprintf(NULL, 0, format, args);// not counting the terminating null character.
        if (msg_len <= 0)return output;

        auto offset = output.length();
        output.resize(offset + (size_t)msg_len);// string[string.size()] is guaranteed to be '\0'
        vsprintf(&output[offset], format, args);
#endif
        if (postfix && *postfix)output.append(postfix);
        return output;
    }

    void set_output_func(void (*func)(const char*))
    {
        print = func;
    }
    void set_exit_func(void(*func)())
    {
        fexit = func;
    }
    void error(const char* format, ...)
    {
#if defined _WIN32
        enable_console_color_for_windows();
#endif

        // format message
        va_list args;
        va_start(args, format);
        std::string msg = format_(format, args, "\n***ERROR*** ", "\n");
        va_end(args);

        // output message to stderr
        std::string msg_with_color = "\u001b[31m" + msg + "\u001b[0m"; // use red text
        fputs(msg_with_color.c_str(), stderr);

        // call user output function
        if (print)print(msg.c_str());

        ASSERT(false);

        if (fexit)fexit(); // call user exit function
        else      throw std::runtime_error(msg);
    }
    void warn(const char* format, ...)
    {
#if defined _WIN32
        enable_console_color_for_windows();
#endif

        // format message
        va_list args;
        va_start(args, format);
        std::string msg = format_(format, args, "\n***WARNING*** ", "\n");
        va_end(args);

        // output message to stdout
        auto msg_with_color = "\u001b[33m" + msg + "\u001b[0m"; // use yellow text
        fputs(msg_with_color.c_str(), stdout);

        // call user output function
        if (print)print(msg.c_str());

        ASSERT(false);
    }
    void info(const char* format, ...)
    {
        // format message
        va_list args;
        va_start(args, format);
        std::string msg = format_(format, args);
        va_end(args);

        // output message to stdout
        fputs(msg.c_str(), stdout);

        // call user output function
        if (print)print(msg.c_str());
    }
    void debug(const char* format, ...)
    {
        // format message
        va_list args;
        va_start(args, format);
        std::string msg = format_(format, args).c_str();
        va_end(args);

        // output message to stdout
        fputs(msg.c_str(), stdout);

        // call user output function
        if (print)print(msg.c_str());
    }
}
