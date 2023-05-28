#pragma once

namespace EasyLib{

    //! @brief Setting the output function for recording message.
    //! @param func Output function.
    void set_output_func(void (*func)(const char*));

    //! @brief Setting the exit function which will be invoked by \error function.
    //! @param func Exit function.
    void set_exit_func(void(*func)());

    //! @brief Writing an error message (automatic adding a prefix of '\n***ERROR***' and postfix of end-of-line).
    //! @param format Pointer to a null-terminated multi-byte string specifying how to interpret the data, see standard \printf function.
    //! @param ...    Arguments specifying data to print.
    //! @note
    //!  + This function will call the exiting function if it is supplied by used through \set_exit_func.
    //!  + Otherwise, a std::runtime_error will be thrown.
    void error(const char* format, ...);

    //! @brief Writing a warning message (automatic adding a prefix of '\n***WARNING***' and postfix of end-of-line).
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
