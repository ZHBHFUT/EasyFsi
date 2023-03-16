//#include <mpi.h>
#include <cstring>

#include "Logger.hpp"
#include "MPICommunicator.hpp"

namespace EasyLib {

    MPICommunicator::MPICommunicator(int mpi_comm, int rank, int size)
        :comm_(mpi_comm), rank_(rank), size_(size)
    {}

    void MPICommunicator::set_constant(const char* name, int value)
    {
        if      (strcmp(name, "MPI_INT16_T"      ) == 0)MPI_INT16_T_ = value;
        else if (strcmp(name, "MPI_INT32_T"      ) == 0)MPI_INT32_T_ = value;
        else if (strcmp(name, "MPI_INT64_T"      ) == 0)MPI_INT64_T_ = value;
        else if (strcmp(name, "MPI_FLOAT"        ) == 0 || strcmp(name, "MPI_REAL4") == 0)MPI_FLOAT_  = value;
        else if (strcmp(name, "MPI_DOUBLE"       ) == 0 || strcmp(name, "MPI_REAL8") == 0)MPI_DOUBLE_ = value;
        else if (strcmp(name, "MPI_CHAR"         ) == 0)MPI_CHAR_ = value;
        else if (strcmp(name, "MPI_SUCCESS"      ) == 0)MPI_SUCCESS_ = value;
        else if (strcmp(name, "MPI_COMM_WORLD"   ) == 0)MPI_COMM_WORLD_ = value;
        else if (strcmp(name, "MPI_COMM_SELF"    ) == 0)MPI_COMM_SELF_ = value;
        else if (strcmp(name, "MPI_COMM_NULL"    ) == 0)MPI_COMM_NULL_ = value;
        else {
            error("unsupported MPI constant \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
        }
    }
    void MPICommunicator::set_constant(const char* name, void* value)
    {
        if (strcmp(name, "MPI_STATUS_IGNORE") == 0)MPI_STATUS_IGNORE_ = (int*)value;
        else {
            error("unsupported MPI constant \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
        }
    }
    void MPICommunicator::set_function(const char* name, void* value)
    {
        if      (strcmp(name, "MPI_Send") == 0)MPI_Send_ = (func_MPI_Send)value;
        else if (strcmp(name, "MPI_Recv") == 0)MPI_Recv_ = (func_MPI_Recv)value;
        else {
            error("unsupported MPI function \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
        }
    }

    bool MPICommunicator::send(const int16_t* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send_(data, count, MPI_INT16_T_, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::send(const int32_t* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send_(data, count, MPI_INT32_T_, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::send(const int64_t* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send_(data, count, MPI_INT64_T_, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::send(const double* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send_(data, count, MPI_DOUBLE_, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::send(const float* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send_(data, count, MPI_FLOAT_, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::send(const char* data, int count, int dest_rank, int tag)
    {
        auto ret = MPI_Send_(data, count, MPI_CHAR_, dest_rank, tag, comm_);
        return ret == MPI_SUCCESS_;
    }

    bool MPICommunicator::recv(int16_t* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv_(data, count, MPI_INT16_T_, src_rank, tag, comm_, MPI_STATUS_IGNORE_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::recv(int32_t* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv_(data, count, MPI_INT32_T_, src_rank, tag, comm_, MPI_STATUS_IGNORE_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::recv(int64_t* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv_(data, count, MPI_INT64_T_, src_rank, tag, comm_, MPI_STATUS_IGNORE_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::recv(double* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv_(data, count, MPI_DOUBLE_, src_rank, tag, comm_, MPI_STATUS_IGNORE_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::recv(float* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv_(data, count, MPI_FLOAT_, src_rank, tag, comm_, MPI_STATUS_IGNORE_);
        return ret == MPI_SUCCESS_;
    }
    bool MPICommunicator::recv(char* data, int count, int src_rank, int tag)
    {
        auto ret = MPI_Recv_(data, count, MPI_CHAR_, src_rank, tag, comm_, MPI_STATUS_IGNORE_);
        return ret == MPI_SUCCESS_;
    }
}
