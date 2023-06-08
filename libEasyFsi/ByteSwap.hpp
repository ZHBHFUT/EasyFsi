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
//! @file       ByteSwap.hpp
//!             The definition of byteswap function.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#ifdef __cpp_lib_byteswap
#include <bit>        // use std::byteswap if is C++23
#else
#include <cstdint>    // uint??_t, int??_t
#include <memory>     // memcpy
#ifdef __linux__
#include <byteswap.h> // __builtin_bswap??
#elif _WIN32
#include <intrin.h>   // _byteswap_??
#else
#error "unsupported platform!"
#endif
#endif

namespace EasyLib {

#ifdef __cpp_lib_byteswap
    template<typename T>
    inline T byteswap(T val) { return std::byteswap(val); }
#else
    inline uint16_t byteswap(uint16_t val)
    {
#ifdef __linux__
        return __builtin_bswap16(val);
#elif _WIN32
        return _byteswap_ushort(val);
#else
        // TBD
#endif
    }

    inline int16_t byteswap(int16_t val)
    {
        uint16_t v;
        std::memcpy(&v, &val, sizeof(val));
        v = byteswap(v);
        std::memcpy(&val, &v, sizeof(val));
        return val;
    }

    inline uint32_t byteswap(uint32_t val)
    {
#ifdef __linux__
        return __builtin_bswap32(val);
#elif _WIN32
        return _byteswap_ulong(val);
#else
        // TBD
#endif
    }

    inline int32_t byteswap(int32_t val)
    {
        uint32_t v;
        std::memcpy(&v, &val, sizeof(val));
        v = byteswap(v);
        std::memcpy(&val, &v, sizeof(val));
        return val;
    }

    inline uint64_t byteswap(uint64_t val)
    {
#ifdef __linux__
        return __builtin_bswap64(val);
#elif _WIN32
        return _byteswap_uint64(val);
#else
        // TBD
#endif
    }

    inline int64_t byteswap(int64_t val)
    {
        uint64_t v;
        std::memcpy(&v, &val, sizeof(val));
        v = byteswap(v);
        std::memcpy(&val, &v, sizeof(val));
        return val;
    }
#endif
}