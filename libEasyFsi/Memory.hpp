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
//! @file       Memory.hpp
//!             The definition memory management functions.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
    typedef void* (*func_allocate)(size_t nbytes);
    typedef void  (*func_deallocate)(void* pointer);

    void set_allocator(func_allocate falloc, func_deallocate fdealloc);

    void* allocate(size_t size_in_bytes);

    void  deallocate(void* pointer);

#ifdef __cplusplus
}
#endif
