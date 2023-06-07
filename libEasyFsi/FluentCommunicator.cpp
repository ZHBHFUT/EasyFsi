#include <cstring>

#include "Logger.hpp"
#include "FluentCommunicator.hpp"

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

}
