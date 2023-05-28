#pragma once
#include "Communicator.hpp"

namespace EasyLib {

    //typedef int(__stdcall *func_MPI_Initialized)(int* flag);
    //typedef int(__stdcall *func_MPI_Init)(const int* argc, char*** argv);
    //typedef int(__stdcall *func_MPI_Comm_disconnect)(int* comm);
    //typedef int(__stdcall *func_MPI_Comm_rank)(int comm, int* rank);
    //typedef int(__stdcall *func_MPI_Comm_size)(int comm, int* size);
    typedef int(__stdcall *func_MPI_Send)(const void* buffer, int count, int datatype, int dest, int tag, int comm);
    typedef int(__stdcall *func_MPI_Recv)(void* buffer, int count, int datatype, int source, int tag, int comm);

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

        void init(int /*argc*/, const char** /*argv*/)final {}

        void disconnect()final {}

        int rank()const noexcept final { return rank_; }
        int size()const noexcept final { return size_; }

        int comm()const noexcept { return comm_; }
        int get_constant(const char* name)const;

        bool send(const int16_t* data, int count, int dest_rank, int tag)final;
        bool send(const int32_t* data, int count, int dest_rank, int tag)final;
        bool send(const int64_t* data, int count, int dest_rank, int tag)final;
        bool send(const double* data, int count, int dest_rank, int tag)final;
        bool send(const float* data, int count, int dest_rank, int tag)final;
        bool send(const char* data, int count, int dest_rank, int tag)final;

        bool recv(int16_t* data, int count, int src_rank, int tag)final;
        bool recv(int32_t* data, int count, int src_rank, int tag)final;
        bool recv(int64_t* data, int count, int src_rank, int tag)final;
        bool recv(double* data, int count, int src_rank, int tag)final;
        bool recv(float* data, int count, int src_rank, int tag)final;
        bool recv(char* data, int count, int src_rank, int tag)final;

    private:
        int comm_{ 0 }, rank_{ 0 }, size_{ 1 };
        func_MPI_Send MPI_Send_{ nullptr };
        func_MPI_Recv MPI_Recv_{ nullptr };
        int  MPI_SUCCESS_      { -1 };
        int  MPI_INT16_T_      { -1 };
        int  MPI_INT32_T_      { -1 };
        int  MPI_INT64_T_      { -1 };
        int  MPI_FLOAT_        { -1 };
        int  MPI_DOUBLE_       { -1 };
        int  MPI_CHAR_         { -1 };
        //int* MPI_STATUS_IGNORE_{ nullptr };
        //int  MPI_COMM_WORLD_   { -1 };
        //int  MPI_COMM_SELF_    { -1 };
        //int  MPI_COMM_NULL_    { -1 };
    };
}
