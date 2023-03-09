#pragma once
#include "Communicator.hpp"

namespace EasyLib {

    //! @brief Fluent communicator, see mport.h
    class FluentCommunicator : public Communicator
    {
    public:
        using func_MPT_csend = void(int mpid, void* data, unsigned int n, int data_type, int tag, const char* file, int line);
        using func_MPT_crecv = int (int mpid, void* data, unsigned int n, int data_type, int tag, const char* file, int line);

        using Communicator::send;
        using Communicator::recv;

        FluentCommunicator() = default;
        FluentCommunicator(const FluentCommunicator&) = default;
        FluentCommunicator& operator = (const FluentCommunicator&) = default;
        virtual ~FluentCommunicator() = default;

        FluentCommunicator(int myid, int np, func_MPT_csend* csend, func_MPT_crecv* crecv);

        void create(int myid, int np, func_MPT_csend* csend, func_MPT_crecv* crecv);

        //! @brief set internal constants, name = one of: MYID,NP,MPT_CHAR,MPT_SHORT,MPT_INT,MPT_LONG,MPT_FLOAT,MPT_DOUBLE,MPT_LONG_LONG
        void set_constant(const char* name, int   value)final;
        void set_constant(const char* name, void* pointer)final;
        //! @brief set internal functions, name = one of: MPT_csend,MPT_crecv
        void set_function(const char* name, void* func_pointer)final;

        inline void set_mpid(int myid, int np) { myid_ = myid; np_ = np; }
        inline void set_MPT_csend(func_MPT_csend* func) { fsend_ = func; }
        inline void set_MPT_crecv(func_MPT_crecv* func) { frecv_ = func; }

        inline void set_MPT_CHAR     (int mpt_char_type = 1) { mpt_char_type_ = mpt_char_type; }
        inline void set_MPT_SHORT    (int mpt_short_type = 2) { mpt_short_type_ = mpt_short_type; }
        inline void set_MPT_INT      (int mpt_int_type = 3) { mpt_int_type_ = mpt_int_type; }
        inline void set_MPT_LONG     (int mpt_long_type = 4) { mpt_long_type_ = mpt_long_type; }
        inline void set_MPT_FLOAT    (int mpt_float_type = 5) { mpt_float_type_ = mpt_float_type; }
        inline void set_MPT_DOUBLE   (int mpt_double_type = 6) { mpt_double_type_ = mpt_double_type; }
        inline void set_MPT_LONG_LONG(int mpt_long_long_type = 8) { mpt_long_long_type_ = mpt_long_long_type; }

        void init(int /*argc*/, const char** /*argv*/)final {}

        void disconnect()final {}

        int rank()const noexcept final { return myid_; }
        int size()const noexcept final { return np_; }

        bool send(const int16_t* data, int count, int dest_rank, int tag)final
        {
            static_assert(sizeof(int16_t) == sizeof(short), "data length not agree");
            fsend_(dest_rank, const_cast<int16_t*>(data), count, mpt_short_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool send(const int32_t* data, int count, int dest_rank, int tag)final
        {
            static_assert(sizeof(int32_t) == sizeof(int), "data length not agree");
            fsend_(dest_rank, const_cast<int32_t*>(data), count, mpt_int_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool send(const int64_t* data, int count, int dest_rank, int tag)final
        {
            static_assert(sizeof(int64_t) == sizeof(long long), "data length not agree");
            fsend_(dest_rank, const_cast<int64_t*>(data), count, mpt_long_long_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool send(const double* data, int count, int dest_rank, int tag)final
        {
            fsend_(dest_rank, const_cast<double*>(data), count, mpt_double_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool send(const float* data, int count, int dest_rank, int tag)final
        {
            fsend_(dest_rank, const_cast<float*>(data), count, mpt_float_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool send(const char* data, int count, int dest_rank, int tag)final
        {
            fsend_(dest_rank, const_cast<char*>(data), count, mpt_char_type_, tag, __FILE__, __LINE__);
            return true;
        }

        bool recv(int16_t* data, int count, int src_rank, int tag)final
        {
            static_assert(sizeof(int16_t) == sizeof(short), "data length not agree");
            frecv_(src_rank, data, count, mpt_short_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool recv(int32_t* data, int count, int src_rank, int tag)final
        {
            static_assert(sizeof(int32_t) == sizeof(int), "data length not agree");
            frecv_(src_rank, data, count, mpt_int_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool recv(int64_t* data, int count, int src_rank, int tag)final
        {
            static_assert(sizeof(int64_t) == sizeof(long long), "data length not agree");
            frecv_(src_rank, data, count, mpt_long_long_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool recv(double* data, int count, int src_rank, int tag)final
        {
            frecv_(src_rank, data, count, mpt_double_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool recv(float* data, int count, int src_rank, int tag)final
        {
            frecv_(src_rank, data, count, mpt_float_type_, tag, __FILE__, __LINE__);
            return true;
        }
        bool recv(char* data, int count, int src_rank, int tag)final
        {
            frecv_(src_rank, data, count, mpt_char_type_, tag, __FILE__, __LINE__);
            return true;
        }

    private:
        int myid_{ -1 };
        int np_  { 0 };

        func_MPT_csend* fsend_{ nullptr };
        func_MPT_crecv* frecv_{ nullptr };

        int mpt_char_type_     { 1 }; // type of int, see MPT_Datatype in mport.h
        int mpt_short_type_    { 2 }; // type of short, see MPT_Datatype in mport.h
        int mpt_int_type_      { 3 }; // type of int, see MPT_Datatype in mport.h
        int mpt_long_type_     { 4 }; // type of long, see MPT_Datatype in mport.h
        int mpt_float_type_    { 5 }; // type of float, see MPT_Datatype in mport.h
        int mpt_double_type_   { 6 }; // type of double, see MPT_Datatype in mport.h
        int mpt_long_long_type_{ 8 }; // type of long long, see MPT_Datatype in mport.h
    };
}
