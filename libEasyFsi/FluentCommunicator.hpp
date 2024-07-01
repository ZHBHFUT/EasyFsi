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
//! @file       FluentCommunicator.hpp
//!             The definition FluentCommunicator class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include "Communicator.hpp"

namespace EasyLib {

    //! @brief Fluent communicator, see mport.h
    class FluentCommunicator : public Communicator
    {
    public:
        using func_MPT_csend = void(int mpid, const void* data, unsigned int n, int data_type, int tag, const char* file, int line);
        using func_MPT_crecv = int (int mpid,       void* data, unsigned int n, int data_type, int tag, const char* file, int line);

        using Communicator::send;
        using Communicator::recv;

        FluentCommunicator() = default;
        FluentCommunicator(const FluentCommunicator&) = default;
        FluentCommunicator& operator = (const FluentCommunicator&) = default;
        virtual ~FluentCommunicator() = default;

        FluentCommunicator(int myid, int np, func_MPT_csend* csend, func_MPT_crecv* crecv);

        //! @brief set internal constants, name = one of: MYID,NP,MPT_CHAR,MPT_SHORT,MPT_INT,MPT_LONG,MPT_FLOAT,MPT_DOUBLE,MPT_LONG_LONG
        void set_constant(const char* name, int   value  )final;
        void set_constant(const char* name, void* pointer)final;

        //! @brief set internal functions, name = one of: MPT_csend,MPT_crecv
        void set_function(const char* name, void* func_pointer)final;

        inline void set_mpid(int myid, int np) { myid_ = myid; np_ = np; }
        inline void set_MPT_csend(func_MPT_csend* func) { fsend_ = func; }
        inline void set_MPT_crecv(func_MPT_crecv* func) { frecv_ = func; }

        inline void set_MPT_CHAR     (int mpt_char_type      = 1) { mpt_char_type_      = mpt_char_type; }
        inline void set_MPT_SHORT    (int mpt_short_type     = 2) { mpt_short_type_     = mpt_short_type; }
        inline void set_MPT_INT      (int mpt_int_type       = 3) { mpt_int_type_       = mpt_int_type; }
        inline void set_MPT_LONG     (int mpt_long_type      = 4) { mpt_long_type_      = mpt_long_type; }
        inline void set_MPT_FLOAT    (int mpt_float_type     = 5) { mpt_float_type_     = mpt_float_type; }
        inline void set_MPT_DOUBLE   (int mpt_double_type    = 6) { mpt_double_type_    = mpt_double_type; }
        inline void set_MPT_LONG_LONG(int mpt_long_long_type = 8) { mpt_long_long_type_ = mpt_long_long_type; }
        inline void set_MPT_UINT     (int mpt_uint_type      = 10) { mpt_uint_type_     = mpt_uint_type; }

        //void init(int /*argc*/, const char** /*argv*/)final {}

        void disconnect()final {}

        int rank()const noexcept final { return myid_; }
        int size()const noexcept final { return np_; }

        void send(const void* data, int count, DataType type, int dest_rank, int tag)final
        {
            fsend_(dest_rank, reinterpret_cast<const int8_t*>(data), count, get_mpt_type(type), tag, __FILE__, __LINE__);
        }

        void recv(void* data, int count, DataType type, int src_rank, int tag)final
        {
            frecv_(src_rank, data, count, get_mpt_type(type), tag, __FILE__, __LINE__);
        }

        void async_send(const void* data, int count, DataType type, int dest_rank, int tag)final;
        void wait()final;

        inline int get_mpt_type(DataType type)const noexcept
        {
            static_assert(sizeof(int8_t ) == sizeof(char), "data length not agree");
            static_assert(sizeof(int16_t) == sizeof(short), "data length not agree");
            static_assert(sizeof(int32_t) == sizeof(int), "data length not agree");
            static_assert(sizeof(int64_t) == sizeof(long long), "data length not agree");

            switch (type) {
            case DataType::cint8:    return mpt_char_type_;
            case DataType::cint16:   return mpt_short_type_;
            case DataType::cint32:   return mpt_int_type_;
            case DataType::cint64:   return mpt_long_long_type_;
            case DataType::cuint8:   return mpt_char_type_;
            case DataType::cuint16:  return mpt_short_type_;
            case DataType::cuint32:  return mpt_uint_type_;
            case DataType::cuint64:  return mpt_long_long_type_;
            case DataType::cdouble:  return mpt_double_type_;
            case DataType::cfloat:   return mpt_float_type_;
            case DataType::cchar:    return mpt_char_type_;
            default:
                return mpt_null_type_;
            }
        }

    private:
        int myid_{ -1 };
        int np_  { 0 };

        func_MPT_csend* fsend_{ nullptr };
        func_MPT_crecv* frecv_{ nullptr };

        int mpt_null_type_     { 0 }; // MPT_DATATYPE_NULL
        int mpt_char_type_     { 1 }; // type of int, see MPT_Datatype in mport.h
        int mpt_short_type_    { 2 }; // type of short, see MPT_Datatype in mport.h
        int mpt_int_type_      { 3 }; // type of int, see MPT_Datatype in mport.h
        int mpt_long_type_     { 4 }; // type of long, see MPT_Datatype in mport.h
        int mpt_float_type_    { 5 }; // type of float, see MPT_Datatype in mport.h
        int mpt_double_type_   { 6 }; // type of double, see MPT_Datatype in mport.h
        int mpt_long_long_type_{ 8 }; // type of long long, see MPT_Datatype in mport.h
        int mpt_uint_type_     { 10 }; // type of unsigned int, see MPT_Datatype in mport.h
    };
}
