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
//! @file       DistributedBoundary.hpp
//!             The definition of DistributedBoundary class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include "Boundary.hpp"
#include "Communicator.hpp"

namespace EasyLib {

    //! @brief Distributed boundary
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
