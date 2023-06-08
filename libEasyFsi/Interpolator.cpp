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
//! @file       Interpolator.cpp
//!             The implement of Interpolator class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <fstream>
#include <sstream>

#include "Field.hpp"
#include "Boundary.hpp"
#include "Logger.hpp"
#include "Interpolator.hpp"

namespace EasyLib {

    template<typename InfoIterator>
    void do_interp_dofs_s2t(int ndof, int_l np_des, InfoIterator it_info, const double** dofs_src, double* dofs_des)
    {
        // interpolate nodal dofs
        for (int_l i = 0; i < np_des; ++i, ++it_info) {
            auto& coeff = *it_info;
            auto  fs    = dofs_src[coeff.src_bd_id]; // source field
            std::fill(dofs_des, dofs_des + ndof, 0);
            for (int j = 0; j < coeff.ndonor; ++j) {
                auto donor  = coeff.donor_nodes  [j];
                auto weight = coeff.donor_weights[j];
                auto dofs_s = fs + ndof * donor;
                for (int k = 0; k < ndof; ++k)
                    dofs_des[k] += weight * dofs_s[k];
            }
            dofs_des += ndof;
        }
    }
    template<typename InfoIterator>
    void do_interp_loads_t2s(int nload, int_l np_des, InfoIterator it_info, const double* load_des, double** load_src)
    {
        // interpolate nodal dofs
        for (int_l i = 0; i < np_des; ++i, ++it_info) {
            auto& coeff = *it_info;
            auto  fs    = load_src[coeff.src_bd_id]; // source loads
            for (int j = 0; j < coeff.ndonor; ++j) {
                auto donor  = coeff.donor_nodes  [j];
                auto weight = coeff.donor_weights[j];
                auto vals_s = fs + nload * donor;
                
                for (int k = 0; k < nload; ++k)
                    vals_s[k] += weight * load_des[k];
            }
            load_des += nload;
        }
    }

    //-----------------------------------------
    // implements of Interpolator
    //-----------------------------------------

    void Interpolator::clear()noexcept
    {
        source_bounds_.clear();
        target_bounds_.clear();
        node_info_.clear();
        face_info_.clear();
        computed_ = false;
    }

    void Interpolator::set_app_id(int source_app, int target_app)noexcept
    {
        source_app_ = source_app;
        target_app_ = target_app;
    }

    int Interpolator::add_source_boundary(Boundary& bd)
    {
        int i = 0;
        for (; i < static_cast<int>(source_bounds_.size()); ++i) {
            if (source_bounds_.at(i) == &bd)return i;
        }
        source_bounds_.push_back(&bd);
        computed_ = false;
        return i;
    }

    int Interpolator::add_target_boundary(Boundary& bd)
    {
        int i = 0;
        for (; i < static_cast<int>(target_bounds_.size()); ++i) {
            if (target_bounds_.at(i) == &bd)return i;
        }
        target_bounds_.push_back(&bd);
        computed_ = false;
        return i;
    }

    void Interpolator::compute_interp_coeff(InterpolationMethod method/* = Automatic*/, int max_donor_for_xps/* = max_donor*/)
    {
        constexpr double eps = 1E-20;

        if (source_bounds_.empty() || target_bounds_.empty())return;

        if (max_donor_for_xps > max_donor) {
            warn("maximum number of donors is great than %d", max_donor);
            max_donor_for_xps = max_donor;
        }
        else if (max_donor_for_xps < min_donor) {
            warn("maximum number of donors is less than %d", min_donor);
            max_donor_for_xps = min_donor;
        }

        // Interpolation method:
        // 
        // If use Automatic:
        //   If face is not defined on source boundary:
        //       If node number <= 1000, on source boundary:
        //           Use GlobalXPS method.
        //       Else:
        //           Use LocalXPS method.
        //   Else:
        //       Use Projection method
        // Else:
        //    Use specified method.
        //

        int_l nn = 0, nf = 0;
        for (auto bs : target_bounds_) {
            nn += bs->node_num();
            nf += bs->face_num();
        }

        node_info_.clear();
        face_info_.clear();
        node_info_.resize(nn);
        face_info_.resize(nf);

        // find nearest source boundary: loop each source bound to accelerating.
        for (int ib = 0; ib < static_cast<int>(source_bounds_.size()); ++ib) {
            auto bd_s = source_bounds_.at(ib);
        
            auto& kdtree = bd_s->kdtree();
        
            auto it_n = node_info_.begin();
            auto it_f = face_info_.begin();
            int_l idx = 0;
            double d2 = 0;

            // loop each target bound and find nearest node
            for (auto bd_t : target_bounds_) {
                for (auto& p : bd_t->node_coords()) {
                    auto& c = *it_n;

                    // skip coincident point
                    if (c.dist_sq <= eps) {
                        ++it_n;
                        continue;
                    }

                    auto n = kdtree.search(p.data(), 1, &idx, &d2);
        
                    // save nearest node
                    if ((n == 1 && c.src_bd_id == invalid_id) || d2 < c.dist_sq) {
                        c.src_bd_id      = ib;
                        c.ndonor         = 1;
                        c.donor_nodes[0] = idx;
                        c.dist_sq        = d2;
                        if (d2 <= eps)c.donor_weights[0] = 1;
                    }
                    // next node
                    ++it_n;
                }
        
                for (auto& p : bd_t->face_centroids()) {
                    auto& c = *it_f;
        
                    // skip coincident point
                    if (c.dist_sq <= eps) {
                        ++it_f;
                        continue;
                    }
        
                    auto n = kdtree.search(p.data(), 1, &idx, &d2);
        
                    // save nearest node
                    if ((n == 1 && it_f->src_bd_id == invalid_id) || d2 < c.dist_sq) {
                        c.src_bd_id      = ib;
                        c.ndonor         = 1;
                        c.donor_nodes[0] = idx;
                        c.dist_sq        = d2;
                        if (d2 <= eps)c.donor_weights[0] = 1;
                    }
                    // next face
                    ++it_f;
                }
            }// next target boundary
        }

        // compute interpolation coefficients
        for (int ib = 0; ib < static_cast<int>(source_bounds_.size()); ++ib) {
            auto bd_s = source_bounds_.at(ib);

            auto it_n = node_info_.begin();
            auto it_f = face_info_.begin();
            
            // use projection method
            if (bd_s->all_high_order() || (bd_s->face_num() > 0 && (method == Projection || method == Automatic))) {
                int_l ids[8] = { 0 };
                double w[8] = { 0 }, d2_sq = 0;
                int   ndonor = 0;

                // loop each target boundary
                for (auto bd_t : target_bounds_) {

                    // loop each node
                    for (auto& p : bd_t->node_coords()) {
                        auto& c = *it_n;

                        // skipping for coincident point.
                        if (c.dist_sq <= eps || c.src_bd_id != ib) { ++it_n; continue; }

                        bd_s->compute_project_interp_coeff(
                            p,
                            ids,
                            w,
                            ndonor,
                            d2_sq   // distance between query point and it's projection.
                        );

                        c.src_bd_id = ib;
                        c.ndonor = ndonor;
                        for (int j = 0; j < ndonor; ++j) {
                            c.donor_nodes[j] = ids[j];
                            c.donor_weights[j] = w[j];
                        }
                        c.dist_sq = d2_sq;

                        ++it_n;
                    }// next node

                    // loop each face-centriod
                    for (auto& p : bd_t->face_centroids()) {
                        auto& c = *it_f;

                        // skipping for coincident point or source bound is not this source.
                        if (c.dist_sq <= eps || c.src_bd_id != ib) { ++it_f; continue; }

                        bd_s->compute_project_interp_coeff(
                            p,
                            ids,
                            w,
                            ndonor,
                            d2_sq  // distance between query point and it's projection.
                        );

                        c.src_bd_id = ib;
                        c.ndonor = ndonor;
                        for (int j = 0; j < ndonor; ++j) {
                            c.donor_nodes[j] = ids[j];
                            c.donor_weights[j] = w[j];
                        }
                        c.dist_sq = d2_sq;

                        ++it_f;
                    } // next face centroid
                } // next target boundary
            }
            // use local spline method
            else {

                // loop each node of target bound
                for (auto bd_t : target_bounds_) {

                    // loop each node of target boundary
                    for (auto& p : bd_t->node_coords()) {
                        auto& c = *it_n;

                        // skipping for coincident point or source bound is not this source.
                        if (c.dist_sq <= eps || c.src_bd_id != ib) { ++it_n; continue; }

                        bd_s->compute_local_xps_interp_coeff(
                            p,
                            max_donor_for_xps,
                            std::span<int_l >(c.donor_nodes, max_donor_for_xps),
                            std::span<double>(c.donor_weights, max_donor_for_xps),
                            c.ndonor
                        );

                        ++it_n;
                    } // next node

                    // loop each face-center of target bound
                    for (auto& p : bd_t->face_centroids()) {
                        auto& c = *it_f;

                        // skipping for coincident point or source bound is not this source.
                        if (c.dist_sq <= eps || c.src_bd_id != ib) { ++it_f; continue; }

                        bd_s->compute_local_xps_interp_coeff(
                            p,
                            max_donor_for_xps,
                            std::span<int_l >(c.donor_nodes, max_donor_for_xps),
                            std::span<double>(c.donor_weights, max_donor_for_xps),
                            c.ndonor
                        );

                        ++it_f;
                    } // next face-centroid
                } // next target boundary
            }
        }

        computed_ = true;
    }

    void Interpolator::interp_dofs_s2t(std::span<Field* const> sources, std::span<Field*> targets)const
    {
        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        // check
        int nerr = 0;
        const int ndof = targets[0]->info->ncomp;
        for (size_t ib = 0; ib < target_bounds_.size(); ++ib) {
            if (targets[ib] == nullptr)continue;

            if (targets[ib]->info->iotype != IncomingDofs) {
                info("***ERROR*** target field \"%s\" is not incoming DOFs!\n", targets[ib]->info->name.c_str());
                ++nerr;
            }
            if (targets[ib]->info->ncomp != ndof) {
                info("***ERROR*** components number of target field \"%s\" not agree!\n", targets[ib]->info->name.c_str());
                ++nerr;
            }
        }
        for (size_t ib = 0; ib < source_bounds_.size(); ++ib) {
            if (sources[ib]->info->iotype != OutgoingDofs) {
                warn("source field \"%s\" is not outgoing DOFs!", sources[ib]->info->name.c_str());
                ++nerr;
            }
            if (sources[ib]->info->location != NodeCentered) {
                warn("source field \"%s\" is not node-centered!", sources[ib]->info->name.c_str());
                ++nerr;
            }
            if (sources[ib]->info->ncomp != ndof) {
                info("***ERROR*** components number of source field \"%s\" not agree!\n", sources[ib]->info->name.c_str());
                ++nerr;
            }
        }
        if (nerr)return;

        std::vector<const double*> sdofs;
        sdofs.reserve(sources.size());
        for (auto f : sources)sdofs.push_back(f->data.data());

        auto nit = node_info_.begin();
        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto&& bt = target_bounds_.at(ib);
            auto nnode = bt->node_num();
            auto nface = bt->face_num();
            auto& tdofs = *targets[ib];

            if (!tdofs.data.empty()) {
                // interpolate nodal dofs
                if (tdofs.info->location == NodeCentered)
                    do_interp_dofs_s2t(ndof, nnode, nit, sdofs.data(), tdofs.data.data());
                // interpolate face dofs
                else
                    do_interp_dofs_s2t(ndof, nface, fit, sdofs.data(), tdofs.data.data());
            }
            
            nit += nnode;
            fit += nface;
        }
    }

    void Interpolator::interp_load_t2s(std::span<Field* const> targets, std::span<Field*> sources/*, bool fill_src_zeros_first = true*/)const
    {
        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        // check
        int nerr = 0;
        const int nload = targets[0]->info->ncomp;
        for (size_t ib = 0; ib < target_bounds_.size(); ++ib) {
            if (targets[ib] == nullptr)continue;

            if (targets[ib]->info->iotype != OutgoingLoads) {
                info("***ERROR*** target field \"%s\" is not outgoing LOADs!\n", targets[ib]->info->name.c_str());
                ++nerr;
            }
            if (targets[ib]->info->ncomp != nload) {
                info("***ERROR*** components number of target field \"%s\" not agree!\n", targets[ib]->info->name.c_str());
                ++nerr;
            }
        }
        for (size_t ib = 0; ib < source_bounds_.size(); ++ib) {
            if (sources[ib]->info->iotype != IncomingLoads) {
                warn("source field \"%s\" is not incoming LOADs!", sources[ib]->info->name.c_str());
                ++nerr;
            }
            if (sources[ib]->info->location != NodeCentered) {
                warn("source field \"%s\" is not node-centered!", sources[ib]->info->name.c_str());
                ++nerr;
            }
            if (sources[ib]->info->ncomp != nload) {
                info("***ERROR*** components number of source field \"%s\" not agree!\n", sources[ib]->info->name.c_str());
                ++nerr;
            }
        }
        if (nerr)return;

        // zero
        //if (fill_src_zeros_first) {
            for (size_t ib = 0; ib < source_bounds_.size(); ++ib) {
                auto& fs = *(sources[ib]);
                fs.data.fill(0);
            }
        //}

        std::vector<double*> sload;
        sload.reserve(sources.size());
        for (auto f : sources)sload.push_back(f->data.data());

        auto nit = node_info_.begin();
        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto&& bt = target_bounds_.at(ib);
            auto nnode = bt->node_num();
            auto nface = bt->face_num();
            auto& tload = *targets[ib];

            if (!tload.data.empty()) {
                // interpolate nodal loads
                if (targets[ib]->info->location == NodeCentered)
                    do_interp_loads_t2s(nload, nnode, nit, tload.data.data(), sload.data());
                // interpolate face loads
                else
                    do_interp_loads_t2s(nload, nface, fit, tload.data.data(), sload.data());
            }

            nit += nnode;
            fit += nface;
        }
    }

    void Interpolator::interp_node_dofs_s2t(int ndof, const double** src_node_dofs, double** des_node_dofs)const
    {
        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        auto nit = node_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto nnode = target_bounds_.at(ib)->node_num();
            
            // interpolate nodal dofs
            if (des_node_dofs[ib])
                do_interp_dofs_s2t(ndof, nnode, nit, src_node_dofs, des_node_dofs[ib]);

            nit += nnode;
        }
    }
    void Interpolator::interp_face_dofs_s2t(int ndof, const double** src_node_dofs, double** des_face_dofs)const
    {
        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto nface = target_bounds_.at(ib)->face_num();
            
            // interpolate face dofs
            if (des_face_dofs[ib])
                do_interp_dofs_s2t(ndof, nface, fit, src_node_dofs, des_face_dofs[ib]);

            fit += nface;
        }
    }
    void Interpolator::interp_node_load_t2s(int nload, double** src_node_load, const double** des_node_load, bool fill_src_zeros_first/* = true*/)const
    {
        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        if (fill_src_zeros_first) {
            for (size_t ib = 0; ib < source_bounds_.size(); ++ib) {
                std::memset(src_node_load[ib], 0, sizeof(double) * nload * source_bounds_.at(ib)->node_num());
            }
        }

        auto nit = node_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto nnode = target_bounds_.at(ib)->node_num();

            if (des_node_load[ib])
                do_interp_loads_t2s(nload, nnode, nit, des_node_load[ib], src_node_load);

            nit += nnode;
        }
    }
    void Interpolator::interp_face_load_t2s(int nload, double** src_node_load, const double** des_face_load, bool fill_src_zeros_first/* = true*/)const
    {
        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        if (fill_src_zeros_first) {
            for (size_t ib = 0; ib < source_bounds_.size(); ++ib) {
                std::memset(src_node_load[ib], 0, sizeof(double) * nload * source_bounds_.at(ib)->node_num());
            }
        }

        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto nface = target_bounds_.at(ib)->face_num();

            if (des_face_load[ib])
                do_interp_loads_t2s(nload, nface, fit, des_face_load[ib], src_node_load);

            fit += nface;
        }
    }

    void Interpolator::interp_all_dofs_s2t()const
    {
        // compute all incoming dofs of target boundaries

        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        std::vector<const double*> sdofs;
        sdofs.reserve(source_bounds_.size());

        auto it_n = node_info_.begin();
        auto it_f = face_info_.begin();
        for (auto&& bt : target_bounds_) {
            auto nnode = bt->node_num();
            auto nface = bt->face_num();

            for (auto&& dofs_t : bt->fields()) {
                if (dofs_t.info->iotype != IncomingDofs)continue;

                sdofs.clear();
                for (auto&& bs : source_bounds_) {
                    auto&& dofs_s = bs->field(dofs_t.info->name.c_str());
                    if      (dofs_s.info->iotype != OutgoingDofs)
                        error("source field is not outgoing DOFs!");
                    else if (dofs_s.info->location != NodeCentered)
                        error("source field is not node-centered!");
                    else if (dofs_s.info->ncomp != dofs_t.info->ncomp)
                        error("component number not agree!");
                    sdofs.push_back(dofs_s.data.data());
                }

                // interpolate
                if (dofs_t.info->location == NodeCentered)
                    do_interp_dofs_s2t(dofs_t.info->ncomp, nnode, it_n, sdofs.data(), const_cast<double*>(dofs_t.data.data()));
                else
                    do_interp_dofs_s2t(dofs_t.info->ncomp, nface, it_f, sdofs.data(), const_cast<double*>(dofs_t.data.data()));
            }
            it_n += nnode;
            it_f += nface;
        }
    }
    void Interpolator::interp_all_load_t2s()const
    {
        // compute all incoming loads of source boundaries
        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        std::vector<double*> sloads;
        sloads.reserve(source_bounds_.size());

        auto nit = node_info_.begin();
        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);
            auto nnode = bt->node_num();
            auto nface = bt->face_num();

            for (auto&& tload : bt->fields()) {
                if (tload.info->iotype != OutgoingLoads)continue;

                sloads.clear();
                for (auto&& bs : source_bounds_) {
                    auto&& sload = bs->field(tload.info->name.c_str());
                    if      (sload.info->iotype != IncomingLoads)
                        error("source field is not incoming load!");
                    else if (sload.info->ncomp != tload.info->ncomp)
                        error("component number not agree!");
                    sloads.push_back(sload.data.data());
                }

                // interpolate nodal loads
                if (tload.info->location == NodeCentered)
                    do_interp_loads_t2s(tload.info->ncomp, nnode, nit, tload.data.data(), sloads.data());
                // interpolate face loads
                else
                    do_interp_loads_t2s(tload.info->ncomp, nface, fit, tload.data.data(), sloads.data());
            }

            nit += nnode;
            fit += nface;
        }
    }

    void Interpolator::interp_dofs_s2t(const char* dof_name)const
    {
        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        std::vector<const double*> sdofs;
        sdofs.reserve(source_bounds_.size());

        auto it_n = node_info_.begin();
        auto it_f = face_info_.begin();
        for (auto&& bt : target_bounds_) {
            auto nnode = bt->node_num();
            auto nface = bt->face_num();

            auto&& dofs_t = bt->field(dof_name);
            if (dofs_t.info->iotype != IncomingDofs)
                error("target field is not incoming DOFs!");

            sdofs.clear();
            for (auto&& bs : source_bounds_) {
                auto&& dofs_s = bs->field(dof_name);
                if      (dofs_s.info->iotype != OutgoingDofs)
                    error("source field is not outgoing DOFs!");
                else if (dofs_s.info->location != NodeCentered)
                    error("source field is not node-centered!");
                else if (dofs_s.info->ncomp != dofs_t.info->ncomp)
                    error("component number not agree!");
                sdofs.push_back(dofs_s.data.data());
            }

            // interpolate
            if (dofs_t.info->location == NodeCentered)
                do_interp_dofs_s2t(dofs_t.info->ncomp, nnode, it_n, sdofs.data(), const_cast<double*>(dofs_t.data.data()));
            else
                do_interp_dofs_s2t(dofs_t.info->ncomp, nface, it_f, sdofs.data(), const_cast<double*>(dofs_t.data.data()));
            
            it_n += nnode;
            it_f += nface;
        }
    }
    void Interpolator::interp_load_t2s(const char* load_name)const
    {
        if (source_bounds_.empty() || target_bounds_.empty())return;
        if (!computed_)error("coefficients are not computed!");

        std::vector<double*> sloads;
        sloads.reserve(source_bounds_.size());

        auto nit = node_info_.begin();
        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);
            auto nnode = bt->node_num();
            auto nface = bt->face_num();

            auto&& tload = bt->field(load_name);
            if (tload.info->iotype != OutgoingLoads)
                error("target field \"%s\" is not outgoing loads!", load_name);

            sloads.clear();
            for (auto&& bs : source_bounds_) {
                auto&& sload = bs->field(load_name);
                if      (sload.info->iotype != IncomingLoads)
                    error("source field is not incoming load!");
                else if (sload.info->ncomp != tload.info->ncomp)
                    error("component number not agree!");
                sloads.push_back(sload.data.data());
            }

            // interpolate nodal loads
            if (tload.info->location == NodeCentered)
                do_interp_loads_t2s(tload.info->ncomp, nnode, nit, tload.data.data(), sloads.data());
            // interpolate face loads
            else
                do_interp_loads_t2s(tload.info->ncomp, nface, fit, tload.data.data(), sloads.data());

            nit += nnode;
            fit += nface;
        }
    }

    void Interpolator::save_coefficients(const char* file)const
    {
        std::ofstream ofs(file, std::ios::binary | std::ios::out);
        if (!ofs.is_open()) { error("failed opening coefficients file: %s", file); return; }
        if (!computed_)return;

        info("\nwriting coefficients file: %s\n", file);

        int buf[2] = {
            static_cast<int>(source_bounds_.size()),
            static_cast<int>(target_bounds_.size())
        };
        ofs.write((const char*)buf, sizeof(buf));

        for (auto bd : source_bounds_) {
            buf[0] = static_cast<int>(bd->node_num());
            buf[1] = static_cast<int>(bd->face_num());
            ofs.write((const char*)buf, sizeof(buf));
        }
        for (auto bd : target_bounds_) {
            buf[0] = static_cast<int>(bd->node_num());
            buf[1] = static_cast<int>(bd->face_num());
            ofs.write((const char*)buf, sizeof(buf));
        }

        if (!node_info_.empty())ofs.write((const char*)node_info_.data(), sizeof(InterpInfo) * node_info_.size());
        if (!face_info_.empty())ofs.write((const char*)face_info_.data(), sizeof(InterpInfo) * face_info_.size());
        ofs.close();

        info("!!!OK!!!\n");
    }

    void Interpolator::load_coefficients(const char* file)
    {
        if (computed_) { warn("coefficient data already exists and will be overrided."); }

        std::ifstream ifs(file, std::ios::binary | std::ios::in);
        if (!ifs.is_open()) { error("failed opening coefficients file: %s", file); return; }

        info("\nreading coefficients file: %s\n", file);

        int buf[2] = { 0 };
        ifs.read((char*)buf, sizeof(buf));
        if (buf[0] != source_bounds_.size() ||
            buf[1] != target_bounds_.size()) {
            error("boundary number not agree");
            return;
        }
        for (auto bd : source_bounds_) {
            ifs.read((char*)buf, sizeof(buf));
            if (buf[0] != bd->node_num() ||
                buf[1] != bd->face_num()) {
                error("node or face number not agree");
                return;
            }
        }
        int nn = 0, nf = 0;
        for (auto bd : target_bounds_) {
            ifs.read((char*)buf, sizeof(buf));
            if (buf[0] != bd->node_num() ||
                buf[1] != bd->face_num()) {
                error("node or face number not agree");
                return;
            }
            nn += buf[0];
            nf += buf[1];
        }

        node_info_.resize(nn);
        face_info_.resize(nf);
        if (nn > 0)ifs.read((char*)node_info_.data(), sizeof(InterpInfo) * node_info_.size());
        if (nf > 0)ifs.read((char*)face_info_.data(), sizeof(InterpInfo) * face_info_.size());

        if (!ifs.good()) { error("failed reading coefficients data"); return; }

        ifs.close();

        computed_ = true;

        info("!!!OK!!!\n");
    }
}
