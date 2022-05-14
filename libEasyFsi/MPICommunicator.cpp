#include <mpi.h>

#include "MPICommunicator.hpp"

namespace EasyLib {

    MPICommunicator::MPICommunicator()
        :comm_(MPI_COMM_NULL)
    {}

    MPICommunicator::MPICommunicator(int mpi_comm)
        :comm_(mpi_comm)
    {}

    void MPICommunicator::create(int mpi_comm)
    {
        comm_ = mpi_comm;
    }

    int MPICommunicator::rank()const noexcept
    {
        int r = 0;
        MPI_Comm_rank(comm_, &r);
        return r;
    }
    int MPICommunicator::size()const noexcept
    {
        int r = 0;
        MPI_Comm_size(comm_, &r);
        return r;
    }

    bool MPICommunicator::send(const int16_t* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send(data, count, MPI_INT16_T, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::send(const int32_t* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send(data, count, MPI_INT32_T, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::send(const int64_t* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send(data, count, MPI_INT64_T, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::send(const double* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send(data, count, MPI_DOUBLE, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::send(const float* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send(data, count, MPI_FLOAT, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::send(const char* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send(data, count, MPI_CHAR, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS;
    }

    bool MPICommunicator::recv(int16_t* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv(data, count, MPI_INT16_T, src_rank, tag, comm_, MPI_STATUS_IGNORE);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::recv(int32_t* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv(data, count, MPI_INT32_T, src_rank, tag, comm_, MPI_STATUS_IGNORE);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::recv(int64_t* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv(data, count, MPI_INT64_T, src_rank, tag, comm_, MPI_STATUS_IGNORE);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::recv(double* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv(data, count, MPI_DOUBLE, src_rank, tag, comm_, MPI_STATUS_IGNORE);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::recv(float* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv(data, count, MPI_FLOAT, src_rank, tag, comm_, MPI_STATUS_IGNORE);
        return ret == MPI_SUCCESS;
    }
    bool MPICommunicator::recv(char* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv(data, count, MPI_CHAR, src_rank, tag, comm_, MPI_STATUS_IGNORE);
        return ret == MPI_SUCCESS;
    }
}
