#pragma once
#include <string>
#include "Index.hpp"

namespace EasyLib {

    class IndexSet;
    class DynamicVector;
    class DynamicMatrix;
    class MeshConnectivity;
    class Boundary;

    class Communicator
    {
    public:
        Communicator() = default;
        virtual ~Communicator() = default;

        virtual int rank()const noexcept = 0;

        virtual int size()const noexcept = 0;

        virtual bool send(const int_l*  data, int count, int dest_rank, int tag) = 0;
        virtual bool send(const int_g*  data, int count, int dest_rank, int tag) = 0;
        virtual bool send(const double* data, int count, int dest_rank, int tag) = 0;
        virtual bool send(const float*  data, int count, int dest_rank, int tag) = 0;
        virtual bool send(const char*   data, int count, int dest_rank, int tag) = 0;

        virtual bool recv(int_l*  data, int count, int src_rank, int tag) = 0;
        virtual bool recv(int_g*  data, int count, int src_rank, int tag) = 0;
        virtual bool recv(double* data, int count, int src_rank, int tag) = 0;
        virtual bool recv(float*  data, int count, int src_rank, int tag) = 0;
        virtual bool recv(char*   data, int count, int src_rank, int tag) = 0;

        virtual void disconnect() = 0;

        bool send(const IndexSet& is, int dest_rank, int tag);

        bool recv(IndexSet& is, int src_rank, int tag);

        bool send(const DynamicVector& dv, int dest_rank, int tag);

        bool recv(DynamicVector& dv, int src_rank, int tag);

        bool send(const DynamicMatrix& dm, int dest_rank, int tag);

        bool recv(DynamicMatrix& dm, int src_rank, int tag);

        bool send(const MeshConnectivity& mc, int dest_rank, int tag);

        bool recv(MeshConnectivity& mc, int src_rank, int tag);

        bool send(const Boundary& bound, int dest_rank, int tag);

        bool recv(Boundary& bound, int src_rank, int tag);

        bool send(const std::string& str, int dest_rank, int tag);

        bool recv(std::string& str, int src_rank, int tag);
    };

}
