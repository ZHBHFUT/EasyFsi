#pragma once

namespace EasyLib{

    //! @brief Setting the output function for recording message.
    //! @param func Output function. The standard \puts function will be used if \func is null.
    void set_output_func(void (*func)(const char*));

    //! @brief Setting the exit function which will be invoked by \error function.
    //! @param func Exit function. The standard \abort function will be used if \func is null.
    void set_exit_func(void(*func)());

    //! @brief Writing an error message (automatic adding a prefix of '***ERROR***' ) and call \abort to exit this program.
    //! @param format Pointer to a null-terminated multi-byte string specifying how to interpret the data, see standard \printf function.
    //! @param ...    Arguments specifying data to print.
    void error(const char* format, ...);

    //! @brief Writing a warning message (automatic adding a prefix of '***WARNING***' ).
    //! @param format Pointer to a null-terminated multi-byte string specifying how to interpret the data, see standard \printf function.
    //! @param ...    Arguments specifying data to print.
    void warn (const char* format, ...);

    //! @brief Writing a common message.
    //! @param format Pointer to a null-terminated multi-byte string specifying how to interpret the data, see standard \printf function.
    //! @param ...    Arguments specifying data to print.
    void info (const char* format, ...);

    //! @brief Writing a debug message.
    //! @param format Pointer to a null-terminated multi-byte string specifying how to interpret the data, see standard \printf function.
    //! @param ...    Arguments specifying data to print.
    void debug(const char* format, ...);
}
