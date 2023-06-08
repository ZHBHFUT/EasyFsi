#pragma once
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
//! @file       Logger.hpp
//!             The definition logging functions.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

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
