#pragma once
#include "Boundary.hpp"
#include "Communicator.hpp"

namespace EasyLib {

    //! @brief 
    class DistributedBoundary
    {
    public:
        using dvec = Boundary::dvec;
        using ivec = Boundary::ivec;
        using vec3 = Boundary::vec3;

        DistributedBoundary() = default;
        DistributedBoundary(const DistributedBoundary&) = default;
        DistributedBoundary& operator = (const DistributedBoundary&) = default;

        void clear();

        void assemble(Boundary& local_bound, Communicator& comm, int root_rank);

        //! @brief Get assembled full boundary data, will be empty if is not root.
        const Boundary& full_boundary()const noexcept { return full_bound_; }
        Boundary& full_boundary() noexcept { return full_bound_; }

        //! @brief Get local boundary.
        const Boundary& local_boundary()const noexcept { return *local_bound_; }

        void gather_node_fields    (int nfields, const double* local_fields, double* global_fields);
        void gather_face_fields    (int nfields, const double* local_fields, double* global_fields);
        void accumulate_node_fields(int nfields, const double* local_fields, double* global_fields);
        void accumulate_face_fields(int nfields, const double* local_fields, double* global_fields);
        void scatter_node_fields   (int nfields, const double* global_fields, double* local_fields);
        void scatter_face_fields   (int nfields, const double* global_fields, double* local_fields);
        
        //void allocate_node_fields  (int nfields, double** local_fields, double** global_fields);
        //void allocate_face_fields  (int nfields, double** local_fields, double** global_fields);
        //void delete_fields         (double** local_fields, double** global_fields);

        //void gather_fields (const char* field_name);
        //void scatter_fields(const char* field_name);

    private:
        Boundary* local_bound_{ nullptr };

        Communicator* comm_{ nullptr };
        int           root_{ 0 };

        // data bellow are available only on root.

        Boundary      full_bound_;
        ivec          part_nodes_ia_g_{ 0 };
        ivec          part_nodes_ja_g_;
        ivec          part_faces_ia_g_{ 0 };

        dvec          buffer_g_;
    };

}
