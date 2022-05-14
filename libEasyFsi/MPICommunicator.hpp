#pragma once
#include "Communicator.hpp"

namespace EasyLib {

    class MPICommunicator : public Communicator
    {
    public:
        using Communicator::send;
        using Communicator::recv;

        MPICommunicator();
        explicit MPICommunicator(int mpi_comm);

        MPICommunicator(const MPICommunicator&) = default;
        MPICommunicator& operator = (const MPICommunicator&) = default;
        virtual ~MPICommunicator() = default;

        void create(int mpi_comm);

        void init(int /*argc*/, const char** /*argv*/)final {}

        int rank()const noexcept override;
        int size()const noexcept override;

        bool send(const int16_t* data, int count, int dest_rank, int tag)override;
        bool send(const int32_t* data, int count, int dest_rank, int tag)override;
        bool send(const int64_t* data, int count, int dest_rank, int tag)override;
        bool send(const double* data, int count, int dest_rank, int tag)override;
        bool send(const float* data, int count, int dest_rank, int tag)override;
        bool send(const char* data, int count, int dest_rank, int tag)override;

        bool recv(int16_t* data, int count, int src_rank, int tag)override;
        bool recv(int32_t* data, int count, int src_rank, int tag)override;
        bool recv(int64_t* data, int count, int src_rank, int tag)override;
        bool recv(double* data, int count, int src_rank, int tag)override;
        bool recv(float* data, int count, int src_rank, int tag)override;
        bool recv(char* data, int count, int src_rank, int tag)override;

    private:
        int comm_;
    };
}
