#include <cstdio>
#include <cstdarg>
#include <string>

#include "Assert.hpp"
#include "Logger.hpp"

namespace EasyLib {

    void default_output(const char* msg)
    {
        printf("%s", msg);
        puts(msg);
    }

    static void (*print)(const char*) = &default_output;
    static void(*fexit)() = nullptr;

    void set_output_func(void (*func)(const char*))
    {
        print = func ? func : &default_output;
    }

    void set_exit_func(void(*func)())
    {
        fexit = func;// ? func : &::abort;
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
        va_list args;
        va_start(args, format);
        print(format_(format, args, "\n***ERROR*** ", "\n").c_str());
        va_end(args);

        ASSERT(false);

        if (fexit)fexit();
        abort();
    }
    void warn(const char* format, ...)
    {
        va_list args;
        va_start(args, format);
        print(format_(format, args, "\n***WARNING*** ", "\n").c_str());
        va_end(args);

        ASSERT(false);
    }
    void info(const char* format, ...)
    {
        va_list args;
        va_start(args, format);
        print(format_(format, args).c_str());
        va_end(args);
    }
    void debug(const char* format, ...)
    {
        va_list args;
        va_start(args, format);
        print(format_(format, args).c_str());
        va_end(args);
    }
}