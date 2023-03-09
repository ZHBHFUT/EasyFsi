//#include <mpi.h>
#include <cstring>

#include "Logger.hpp"
#include "MPICommunicator.hpp"

namespace EasyLib {

    MPICommunicator::MPICommunicator()
        :comm_(-1)
    {}

    MPICommunicator::MPICommunicator(int mpi_comm)
        :comm_(mpi_comm)
    {}

    void MPICommunicator::create(int mpi_comm)
    {
        comm_ = mpi_comm;

        //if (!MPI_Initialized_    )error("function MPI_Initialized is not specified!");
        //if (!MPI_Init_           )error("function MPI_Init is not specified!");
        //if (!MPI_Comm_disconnect_)error("function MPI_Comm_disconnect is not specified!");
        //if (!MPI_Comm_rank_      )error("function MPI_Comm_rank is not specified!");
        //if (!MPI_Comm_size_      )error("function MPI_Comm_size is not specified!");
        //if (!MPI_Send_           )error("function MPI_Send is not specified!");
        //if (!MPI_Recv_           )error("function MPI_Recv is not specified!");
        //
        //if (MPI_INT16_T_ < 0)error("constant MPI_INT16_T is not specified!");
        //if (MPI_INT32_T_ < 0)error("constant MPI_INT32_T is not specified!");
        //if (MPI_FLOAT_   < 0)error("constant MPI_FLOAT is not specified!");
        //if (MPI_DOUBLE_  < 0)error("constant MPI_DOUBLE is not specified!");
        //if (MPI_CHAR_    < 0)error("constant MPI_CHAR is not specified!");
        //
        //if (!MPI_STATUS_IGNORE_)error("constant MPI_STATUS_IGNORE is not specified!");
        //if (MPI_COMM_WORLD_    < 0)error("constant MPI_COMM_WORLD is not specified!");
        //if (MPI_COMM_SELF_     < 0)error("constant MPI_COMM_SELF is not specified!");
        //if (MPI_COMM_NULL_     < 0)error("constant MPI_COMM_NULL is not specified!");
    }

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
        if      (strcmp(name, "MPI_Initialized"    ) == 0)MPI_Initialized_     = (func_MPI_Initialized)value;
        else if (strcmp(name, "MPI_Init"           ) == 0)MPI_Init_            = (func_MPI_Init)value;
        else if (strcmp(name, "MPI_Comm_disconnect") == 0)MPI_Comm_disconnect_ = (func_MPI_Comm_disconnect)value;
        else if (strcmp(name, "MPI_Comm_rank"      ) == 0)MPI_Comm_rank_       = (func_MPI_Comm_rank)value;
        else if (strcmp(name, "MPI_Comm_size"      ) == 0)MPI_Comm_size_       = (func_MPI_Comm_size)value;
        else if (strcmp(name, "MPI_Send"           ) == 0)MPI_Send_            = (func_MPI_Send)value;
        else if (strcmp(name, "MPI_Recv"           ) == 0)MPI_Recv_            = (func_MPI_Recv)value;
        else {
            error("unsupported MPI function \"%s\", FILE = %s, FUNC = %s, LINE = %s.", name, __FILE__, __FUNCTION__, __LINE__);
        }
    }

    void MPICommunicator::init(int argc, const char** argv)
    {
        int flag = 0;
        MPI_Initialized_(&flag);
        if (!flag)MPI_Init_(&argc, const_cast<char***>(&argv));
    }
    void MPICommunicator::disconnect()
    {
        if (comm_ != MPI_COMM_WORLD_ &&
            comm_ != MPI_COMM_SELF_ &&
            comm_ != MPI_COMM_NULL_)
            MPI_Comm_disconnect_(&comm_);
        else
            comm_ = MPI_COMM_NULL_;
    }
    int MPICommunicator::rank()const noexcept
    {
        int r = -1;
        MPI_Comm_rank_(comm_, &r);
        return r;
    }
    int MPICommunicator::size()const noexcept
    {
        int r = 0;
        MPI_Comm_size_(comm_, &r);
        return r;
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
