#include <cstdio>
#include <cstdarg>
#include <string>

#include <stdexcept>

#include "Assert.hpp"
#include "Logger.hpp"

extern void enable_console_color();

namespace EasyLib {
    static void (*print)(const char*) = nullptr;
    static void (*fexit)() = nullptr;

    void set_output_func(void (*func)(const char*))
    {
        print = func;
    }

    void set_exit_func(void(*func)())
    {
        fexit = func;
    }

    std::string format_(const char* format, va_list args, const char* prefix = nullptr, const char* postfix = nullptr)
    {
        std::string output((prefix && *prefix) ? prefix : "");
        
#ifdef _MSC_VER
        //use _vscprintf, vsprintf_s
        auto msg_len = _vscprintf(format, args);// not counting the terminating null character.
        if (msg_len <= 0)return output;

        auto offset = output.length();
        output.resize(offset + (size_t)msg_len);// string[string.size()] is guaranteed to be '\0'
        //vsprintf(output.data(), fmt, arg_ptr);//C4996
        vsprintf_s(output.data() + offset, msg_len + 1, format, args);
#else
        //use vsnprintf, vsprintf
        auto msg_len = vsnprintf(NULL, 0, format, args);// not counting the terminating null character.
        if (msg_len <= 0)return output;

        auto offset = output.length();
        output.resize(offset + (size_t)msg_len);// string[string.size()] is guaranteed to be '\0'
        vsprintf(output.data() + offset, format, args);
#endif
        if (postfix && *postfix)output.append(postfix);
        return output;
    }

    void error(const char* format, ...)
    {
        enable_console_color();

        va_list args;
        va_start(args, format);
        std::string msg = format_(format, args, "\n***ERROR*** ", "\n");
        va_end(args);

        std::string msg_with_color = "\u001b[31m" + msg + "\u001b[0m"; // use red text
        fputs(msg_with_color.c_str(), stderr);
        if (print)print(msg.c_str());

        ASSERT(false);

        if (fexit)
            fexit();
        else
            throw std::runtime_error(msg);
    }
    void warn(const char* format, ...)
    {
        enable_console_color();

        va_list args;
        va_start(args, format);
        std::string msg = format_(format, args, "\n***WARNING*** ", "\n");
        va_end(args);

        auto msg_with_color = "\u001b[33m" + msg + "\u001b[0m"; // use yellow text
        fputs(msg_with_color.c_str(), stdout);
        if (print)print(msg.c_str());

        ASSERT(false);
    }
    void info(const char* format, ...)
    {
        va_list args;
        va_start(args, format);
        std::string msg = format_(format, args);
        va_end(args);

        fputs(msg.c_str(), stdout);
        if (print)print(msg.c_str());
    }
    void debug(const char* format, ...)
    {
        va_list args;
        va_start(args, format);
        std::string msg = format_(format, args).c_str();
        va_end(args);

        fputs(msg.c_str(), stdout);
        if (print)print(msg.c_str());
    }
}
