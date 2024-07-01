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
//! @file       FluentCommunicator.cpp
//!             The implement of FluentCommunicator class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------
#include "FluentCommunicator.hpp"
#include <cstring>
#include "Logger.hpp"

namespace EasyLib {

    FluentCommunicator::FluentCommunicator(int myid, int np, func_MPT_csend* csend, func_MPT_crecv* crecv)
        :myid_(myid), np_(np), fsend_(csend), frecv_(crecv)
    {}

    void FluentCommunicator::set_constant(const char* name, int value)
    {
        if      (_stricmp(name, "myid"         ) == 0 || _stricmp(name, "rank") == 0)myid_ = value;
        else if (_stricmp(name, "np"           ) == 0 || _stricmp(name, "size") == 0)np_   = value;
        else if (_stricmp(name, "MPT_CHAR"     ) == 0)mpt_char_type_      = value;
        else if (_stricmp(name, "MPT_SHORT"    ) == 0)mpt_short_type_     = value;
        else if (_stricmp(name, "MPT_INT"      ) == 0)mpt_int_type_       = value;
        else if (_stricmp(name, "MPT_LONG"     ) == 0)mpt_long_type_      = value;
        else if (_stricmp(name, "MPT_FLOAT"    ) == 0)mpt_float_type_     = value;
        else if (_stricmp(name, "MPT_DOUBLE"   ) == 0)mpt_double_type_    = value;
        else if (_stricmp(name, "MPT_LONG_LONG") == 0)mpt_long_long_type_ = value;
        else if (_stricmp(name, "MPT_UNSIGNED_INT") == 0)mpt_uint_type_ = value;
        else {
            error("unsupported constant \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
        }
    }
    void FluentCommunicator::set_constant(const char* /*name*/, void* /*pointer*/)
    {
        // do nothing
    }
    void FluentCommunicator::set_function(const char* name, void* func_pointer)
    {
        if      (strcmp(name, "MPT_csend") == 0)fsend_ = (func_MPT_csend*)func_pointer;
        else if (strcmp(name, "MPT_crecv") == 0)frecv_ = (func_MPT_crecv*)func_pointer;
        else {
            error("unsupported function \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
        }
    }

    void FluentCommunicator::async_send(const void* /*data*/, int /*count*/, DataType /*type*/, int /*dest_rank*/, int /*tag*/)
    {
        error("this feature is not implemented! FILE = %s, FUNC = %s, LINE = %s.", __FILE__, __FUNCTION__, __LINE__);
    }
    void FluentCommunicator::wait()
    {
        error("this feature is not implemented! FILE = %s, FUNC = %s, LINE = %s.", __FILE__, __FUNCTION__, __LINE__);
    }
}
