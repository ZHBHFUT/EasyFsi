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
//! @file       MPICommunicator.cpp
//!             The implement of MPICommunicator class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <cstring>

#include "Logger.hpp"
#include "MPICommunicator.hpp"

namespace EasyLib {

    MPICommunicator::MPICommunicator(int mpi_comm, int rank, int size)
        :comm_(mpi_comm), rank_(rank), size_(size)
    {}

    void MPICommunicator::set_constant(const char* name, int value)
    {
        if      (strcmp(name, "MPI_DATATYPE_NULL") == 0)MPI_NULL_T_  = value;
        else if (strcmp(name, "MPI_INT16_T"      ) == 0)MPI_INT16_T_ = value;
        else if (strcmp(name, "MPI_INT32_T"      ) == 0)MPI_INT32_T_ = value;
        else if (strcmp(name, "MPI_INT64_T"      ) == 0)MPI_INT64_T_ = value;
        else if (strcmp(name, "MPI_FLOAT"        ) == 0 || strcmp(name, "MPI_REAL4") == 0)MPI_FLOAT_  = value;
        else if (strcmp(name, "MPI_DOUBLE"       ) == 0 || strcmp(name, "MPI_REAL8") == 0)MPI_DOUBLE_ = value;
        else if (strcmp(name, "MPI_CHAR"         ) == 0)MPI_CHAR_ = value;
        //else if (strcmp(name, "MPI_SUCCESS"      ) == 0)MPI_SUCCESS_ = value;
        //else if (strcmp(name, "MPI_COMM_WORLD"   ) == 0)MPI_COMM_WORLD_ = value;
        //else if (strcmp(name, "MPI_COMM_SELF"    ) == 0)MPI_COMM_SELF_ = value;
        //else if (strcmp(name, "MPI_COMM_NULL"    ) == 0)MPI_COMM_NULL_ = value;
        else {
            error("unsupported MPI constant \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
        }
    }
    //void MPICommunicator::set_constant(const char* name, void* value)
    //{
    //    if (strcmp(name, "MPI_STATUS_IGNORE") == 0)MPI_STATUS_IGNORE_ = (int*)value;
    //    else {
    //        error("unsupported MPI constant \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
    //    }
    //}
    void MPICommunicator::set_function(const char* name, void* value)
    {
        if      (strcmp(name, "MPI_Send") == 0)MPI_Send_ = (func_MPI_Send)value;
        else if (strcmp(name, "MPI_Recv") == 0)MPI_Recv_ = (func_MPI_Recv)value;
        else {
            error("unsupported MPI function \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
        }
    }

    int MPICommunicator::get_constant(const char* name)const
    {
        if      (strcmp(name, "MPI_DATATYPE_NULL") == 0)return MPI_NULL_T_;
        else if (strcmp(name, "MPI_INT16_T"   ) == 0)return MPI_INT16_T_;
        else if (strcmp(name, "MPI_INT32_T"   ) == 0)return MPI_INT32_T_;
        else if (strcmp(name, "MPI_INT64_T"   ) == 0)return MPI_INT64_T_;
        else if (strcmp(name, "MPI_FLOAT"     ) == 0 || strcmp(name, "MPI_REAL4") == 0)return MPI_FLOAT_;
        else if (strcmp(name, "MPI_DOUBLE"    ) == 0 || strcmp(name, "MPI_REAL8") == 0)return MPI_DOUBLE_;
        else if (strcmp(name, "MPI_CHAR"      ) == 0)return MPI_CHAR_;
        //else if (strcmp(name, "MPI_SUCCESS"   ) == 0)return MPI_SUCCESS_;
        //else if (strcmp(name, "MPI_COMM_WORLD") == 0)return MPI_COMM_WORLD_;
        //else if (strcmp(name, "MPI_COMM_SELF" ) == 0)return MPI_COMM_SELF_;
        //else if (strcmp(name, "MPI_COMM_NULL" ) == 0)return MPI_COMM_NULL_;
        else {
            error("unknown MPI constant \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
            return 0;
        }
    }

    int MPICommunicator::dt2mpt_(DataType type)const noexcept
    {
        static_assert(sizeof(int8_t) == sizeof(char));

        switch (type) {
        case DataType::cint8: [[fallthrough]];
        case DataType::cuint8: [[fallthrough]];
        case DataType::cchar: [[fallthrough]];
        case DataType::cuchar:return MPI_CHAR_;
        case DataType::cint16: [[fallthrough]];
        case DataType::cuint16:return MPI_INT16_T_;
        case DataType::cint32: [[fallthrough]];
        case DataType::cuint32:return MPI_INT32_T_;
        case DataType::cint64: [[fallthrough]];
        case DataType::cuint64:return MPI_INT64_T_;
        case DataType::cfloat:return MPI_FLOAT_;
        case DataType::cdouble:return MPI_DOUBLE_;
        default:
            return MPI_NULL_T_;
        }
    }

    void MPICommunicator::send(const void* data, int count, DataType type, int dest_rank, int tag)
    {
        MPI_Send_(data, count, dt2mpt_(type), dest_rank, tag, comm_);
    }
    void MPICommunicator::recv(void* data, int count, DataType type, int src_rank, int tag)
    {
        MPI_Recv_(data, count, dt2mpt_(type), src_rank, tag, comm_);
    }
}
