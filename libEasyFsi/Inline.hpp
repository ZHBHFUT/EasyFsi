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
//! @file       Inline.hpp
//!             The definition _force_inline_ macro.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#if defined(__clang__) || defined(__GNUC__)
#define _force_inline_ __attribute__((always_inline)) inline
#define _force_inline_ [[gnu::always_inline]] inline
#elif defined(_MSC_VER)
#ifdef _DEBUG
#define _force_inline_ inline
#else
#pragma warning(error: 4714)
#define _force_inline_ __forceinline
#endif
#else
#error Unsupported compiler
#endif
