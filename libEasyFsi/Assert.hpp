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
//! @file       Assert.hpp
//!             The definition of ASSERT macro.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2019, all rights reserved.
//! @data       2019-05-13
//!-------------------------------------------------------------

#ifndef ASSERT
#  ifndef __CUDACC__
#    ifdef _MSC_VER // for vc++
#      include <crtdbg.h> // _ASSERT
#      ifdef _DEBUG
#        define ASSERT(expr) _ASSERT(expr)
#      else
#        define ASSERT(expr)
#      endif
#      define BREAK_IF(expr) (void) (!!(expr)) || (__debugbreak(), 0))
#    else // for other compiler
#      ifdef NDEBUG
#        define ASSERT(expr)
#      else
#        include <assert.h>
#        define ASSERT(expr) assert(expr)
#      endif
#      define BREAK_IF(expr)
#    endif // _MSC_VER
#  else
#    define ASSERT(expr)
#    define BREAK_IF(expr)
#  endif // __CUDACC__
#endif
