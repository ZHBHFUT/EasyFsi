#pragma once
#include "Communicator.hpp"

namespace EasyLib {

    typedef int(__stdcall *func_MPI_Send)(const void* buffer, int count, int datatype, int dest  , int tag, int comm);
    typedef int(__stdcall *func_MPI_Recv)(      void* buffer, int count, int datatype, int source, int tag, int comm);

    class MPICommunicator : public Communicator
    {
    public:
        using Communicator::send;
        using Communicator::recv;

        MPICommunicator() = default;
        explicit MPICommunicator(int mpi_comm, int rank, int size);

        MPICommunicator(const MPICommunicator&) = default;
        MPICommunicator& operator = (const MPICommunicator&) = default;
        virtual ~MPICommunicator() = default;

        //! @brief  set internal constants, name = one of:
        //!   MPI_INT16_T,MPI_INT32_T,MPI_INT64_T,MPI_FLOAT,MPI_DOUBLE,MPI_CHAR
        //!   MPI_SUCCESS
        void set_constant(const char* name, int   value)final;

        //! @brief set internal pointer, name = MPI_STATUS_IGNORE
        void set_constant(const char* /*name*/, void* /*pointer*/)final {}

        //! @brief set internal function pointer, name = one of
        //!   MPI_Send
        //!   MPI_Recv
        void set_function(const char* name, void* func_pointer)final;

        //void init(int /*argc*/, const char** /*argv*/)final {}

        void disconnect()final {}

        int rank()const noexcept final { return rank_; }
        int size()const noexcept final { return size_; }

        int comm()const noexcept { return comm_; }
        int get_constant(const char* name)const;

        void send(const void* data, int count, DataType type, int dest_rank, int tag)final;
        void recv(      void* data, int count, DataType type, int src_rank,  int tag)final;

    private:
        int dt2mpt_(DataType type)const noexcept;

    private:
        int comm_{ 0 }, rank_{ 0 }, size_{ 1 };
        func_MPI_Send MPI_Send_{ nullptr };
        func_MPI_Recv MPI_Recv_{ nullptr };
        //int  MPI_SUCCESS_      { -1 };
        int  MPI_NULL_T_       { -1 };
        int  MPI_INT16_T_      { -1 };
        int  MPI_INT32_T_      { -1 };
        int  MPI_INT64_T_      { -1 };
        int  MPI_FLOAT_        { -1 };
        int  MPI_DOUBLE_       { -1 };
        int  MPI_CHAR_         { -1 };
    };
}
