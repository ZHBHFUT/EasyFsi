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
//! @file       Interpolator.hpp
//!             The definition Interpolator class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <limits>

#include "DynamicMatrix.hpp"

namespace EasyLib {

    class Boundary;
    struct Field;

    enum InterpolationMethod
    {
        //GlobalXPS,  //! Use global IPS/TPS/Spline method.
        LocalXPS,   //! Use local IPS/TPS/Spline method.
        Projection, //! Used projection method, faces should be defined on boundary.
        Mapping,    //! This is not implemented now.
        Automatic   //! Automatic select above method.
    };

    class Interpolator
    {
    public:
        inline static constexpr const int max_donor = 20;
        inline static constexpr const int min_donor = 8;

        struct InterpInfo
        {
            int    src_bd_id{ 0 };
            int    ndonor{ 0 };
            int_l  donor_nodes  [max_donor]{ 0 };
            double donor_weights[max_donor]{ 0 };
            size_t donor_beg{ 0 }, donor_end{ 0 };
            double dist_sq{ std::numeric_limits<double>::max() };
        };

        Interpolator() = default;

        virtual ~Interpolator() = default;

        void clear()noexcept;

        void set_app_id(int source_app, int target_app)noexcept;

        //! @brief Add a source boundary to interpolator. All fields of source boundary must be node-centered.
        //! @param bd The source boundary
        //! @return index in the interpolator.
        int add_source_boundary(Boundary& bd);

        //! @brief Add a target boundary to interpolator.
        //! @param bd The target boundary
        //! @return index in the interpolator.
        int add_target_boundary(Boundary& bd);

        //! @brief Compute interpolation coefficients.
        void compute_interp_coeff(InterpolationMethod method = Automatic, int max_donor_for_xps = max_donor);

        //! @brief Compute incoming DOFs from source bounds to target bounds.
        //! @param sources Field of source boundaries. size = source bounds number.
        //! @param targets Field of target boundaries. size = target bounds number.
        //! @note each of \sources field must be node-centered.
        void interp_dofs_s2t(Span<Field* const> sources, Span<Field*> targets)const;

        //! @brief Compute incoming LOADs from target bounds to source bounds. 
        //! @param targets Field of target boundaries. size = target bounds number.
        //! @param sources Field of source boundaries. size = source bounds number.
        //! @param fill_src_zeros_first  Whether or not to fill the source fields with zeros before computing.
        //! @note
        //!  + each of \sources field must be node-centered.
        void interp_load_t2s(Span<Field* const> targets, Span<Field*> sources/*, bool fill_src_zeros_first = true*/)const;

        //! @brief Compute nodal DOFs of target boundaries by source boundaries.
        //! @param [in]  ndof           DOF number
        //! @param [in]  src_node_dofs  Nodal DOFs of source boundaries: {{u1,v1,...,u2,v2,...}, {u1,v1,...,u2,v2,...},...}
        //! @param [out] des_node_dofs  Nodal DOFs of target boundaries: {{u1,v1,...,u2,v2,...}, {u1,v1,...,u2,v2,...},...}
        void interp_node_dofs_s2t(int ndof, const double** src_node_dofs, double** des_node_dofs)const;
        //! @brief Compute face-centered DOFs of target boundaries by nodal DOFs of source boundaries.
        //! @param [in]  ndof           DOF number
        //! @param [in]  src_node_dofs  Nodal DOFs of source boundaries: {{u1,v1,...,u2,v2,...}, {u1,v1,...,u2,v2,...},...}
        //! @param [out] des_node_dofs  Nodal DOFs of target boundaries: {{u1,v1,...,u2,v2,...}, {u1,v1,...,u2,v2,...},...}
        void interp_face_dofs_s2t(int ndof, const double** src_node_dofs, double** des_face_dofs)const;

        //! @brief Compute nodal loads of source boundaries by target boundaries.
        //! @param [in]  nload          Load number
        //! @param [out] src_node_load  Nodal loads of source boundaries: {{f1,f2,...,g1,g2,...},{f1,f2,...,g1,g2,...},...}
        //! @param [in]  des_node_load  Nodal loads of target boundaries: {{f1,f2,...,g1,g2,...},{f1,f2,...,g1,g2,...},...}
        //! @param [in]  fill_src_zeros_first Should we fill output array with zeros?
        void interp_node_load_t2s(int nload, double** src_node_load, const double** des_node_load, bool fill_src_zeros_first = true)const;

        //! @brief Compute face-centered loads of source boundaries by nodal loads of target boundaries.
        //! @param [in]  nload          Load number
        //! @param [out] src_node_load  Face-centered loads of source boundaries: {{f1,f2,...,g1,g2,...},{f1,f2,...,g1,g2,...},...}
        //! @param [in]  des_face_load  Face-centered loads of target boundaries: {{f1,f2,...,g1,g2,...},{f1,f2,...,g1,g2,...},...}
        //! @param [in]  fill_src_zeros_first Should we fill output array with zeros?
        void interp_face_load_t2s(int nload, double** src_node_load, const double** des_face_load, bool fill_src_zeros_first = true)const;

        void interp_all_dofs_s2t()const;
        void interp_all_load_t2s()const;

        void interp_dofs_s2t(const char* dof_name)const;
        void interp_load_t2s(const char* load_name)const;

        inline int source_app_id()const noexcept { return source_app_; }
        inline int target_app_id()const noexcept { return target_app_; }

        void save_coefficients(const char* file)const;
        void load_coefficients(const char* file);

    private:
        int source_app_{ -1 };
        int target_app_{ -1 };
        std::vector<Boundary*> source_bounds_;
        std::vector<Boundary*> target_bounds_;
        bool computed_{ false };

        std::vector<InterpInfo> node_info_;
        std::vector<InterpInfo> face_info_;
    };
}
