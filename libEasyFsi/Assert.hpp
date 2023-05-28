#pragma once
//!-------------------------------------------------------------
//! @file       Assert.hpp
//!             The defination of ASSERT macro for debugging.
//! @author     ZHANGBing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2019, all rights reserved.
//! @data       2019-05-13
//! @note       Only used for HFUT-CFD group.
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
