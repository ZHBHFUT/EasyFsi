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
//! @file       Communicator.hpp
//!             The definition of Communicator class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <cstdint>
#include <string>

#include "Inline.hpp"
#include "Assert.hpp"

namespace EasyLib {

    class IndexSet;
    class DynamicVector;
    class DynamicMatrix;
    class MeshConnectivity;
    class Boundary;
    struct FieldInfo;

    enum class DataType
    {
        cint8,
        cint16,
        cint32,
        cint64,
        cuint8,
        cuint16,
        cuint32,
        cuint64,
        cfloat,
        cdouble,
        cchar,
        cuchar
    };

    template<typename T>
    struct DataTypeTraits;
    //template<typename T>
    //struct DataTypeTraits
    //{
    //    static_const bool is_basic = false;
    //    static_const int  nbytes   = sizeof(T);
    //};
    template<> struct DataTypeTraits<char>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(char);
        static_const DataType dtype    = DataType::cchar;
    };
    template<> struct DataTypeTraits<int8_t>
    {
        static_const bool     is_basic  = true;
        static_const int      nbytes    = sizeof(int8_t);
        static_const DataType dtype     = DataType::cint8;
    };
    template<> struct DataTypeTraits<int16_t>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(int16_t);
        static_const DataType dtype    = DataType::cint16;
    };
    template<> struct DataTypeTraits<int32_t>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(int32_t);
        static_const DataType dtype    = DataType::cint32;
    };
    template<> struct DataTypeTraits<int64_t>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(int64_t);
        static_const DataType dtype    = DataType::cint64;
    };
    template<> struct DataTypeTraits<uint8_t>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(uint8_t);
        static_const DataType dtype    = DataType::cuint8;
    };
    template<> struct DataTypeTraits<uint16_t>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(uint16_t);
        static_const DataType dtype    = DataType::cuint16;
    };
    template<> struct DataTypeTraits<uint32_t>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(uint32_t);
        static_const DataType dtype    = DataType::cuint32;
    };
    template<> struct DataTypeTraits<uint64_t>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(uint64_t);
        static_const DataType dtype    = DataType::cuint64;
    };
    template<> struct DataTypeTraits<float>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(float);
        static_const DataType dtype    = DataType::cfloat;
    };
    template<> struct DataTypeTraits<double>
    {
        static_const bool     is_basic = true;
        static_const int      nbytes   = sizeof(double);
        static_const DataType dtype    = DataType::cdouble;
    };
    
    inline constexpr int nbyte_of_type(DataType type)
    {
        switch (type) {
        case DataType::cint8:   return sizeof(int8_t);
        case DataType::cint16:  return sizeof(int16_t);
        case DataType::cint32:  return sizeof(int32_t);
        case DataType::cint64:  return sizeof(int64_t);
        case DataType::cuint8:  return sizeof(uint8_t);
        case DataType::cuint16: return sizeof(uint16_t);
        case DataType::cuint32: return sizeof(uint32_t);
        case DataType::cuint64: return sizeof(uint64_t);
        case DataType::cfloat:  return sizeof(float);
        case DataType::cdouble: return sizeof(double);
        case DataType::cchar:   return sizeof(char);
        case DataType::cuchar:  return sizeof(unsigned char);
        default: ASSERT(false); return 0;
        }
    }

    class Communicator
    {
    public:
        
        Communicator() = default;
        virtual ~Communicator() = default;

        virtual void set_constant(const char* name, int   value  ) = 0;
        virtual void set_constant(const char* name, void* pointer) = 0;
        virtual void set_function(const char* name, void* func   ) = 0;

        virtual int rank()const noexcept = 0;
        virtual int size()const noexcept = 0;

        virtual void send(const void* data, int count, DataType type, int dest_rank, int tag) = 0;
        virtual void recv(      void* data, int count, DataType type, int src_rank,  int tag) = 0;

        virtual void async_send(const void* data, int count, DataType type, int dest_rank, int tag) = 0;
        virtual void wait() = 0;

        virtual void disconnect() = 0;

        template<typename T>
        void send(const T* data, int count, int dest_rank, int tag)
        {
            static_assert(DataTypeTraits<T>::is_basic, "invalid data type!");
            send(data, count, DataTypeTraits<T>::dtype, dest_rank, tag);
        }
        template<typename T>
        void recv(T* data, int count, int src_rank, int tag)
        {
            static_assert(DataTypeTraits<T>::is_basic, "invalid data type!");
            recv(data, count, DataTypeTraits<T>::dtype, src_rank, tag);
        }

        void send(const IndexSet& is, int dest_rank, int tag);
        void recv(IndexSet& is, int src_rank, int tag);
        void send(const DynamicVector& dv, int dest_rank, int tag);
        void recv(DynamicVector& dv, int src_rank, int tag);
        void send(const DynamicMatrix& dm, int dest_rank, int tag);
        void recv(DynamicMatrix& dm, int src_rank, int tag);
        void send(const MeshConnectivity& mc, int dest_rank, int tag);
        void recv(MeshConnectivity& mc, int src_rank, int tag);
        void send(const Boundary& bound, int dest_rank, int tag);
        void recv(Boundary& bound, int src_rank, int tag);
        void send(const std::string& str, int dest_rank, int tag);
        void recv(std::string& str, int src_rank, int tag);
        void send(const FieldInfo& info, int dest_rank, int tag);
        void recv(FieldInfo& info, int src_rank, int tag);
    };

}
