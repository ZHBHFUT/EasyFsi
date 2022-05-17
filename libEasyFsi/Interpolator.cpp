#include "Field.hpp"
#include "Boundary.hpp"
#include "Logger.hpp"
#include "Interpolator.hpp"

namespace EasyLib {

    //-----------------------------------------
    // implements of Interpolator
    //-----------------------------------------

    void Interpolator::clear()
    {
        source_bounds_.clear();
        target_bounds_.clear();
        computed_ = false;
    }

    void Interpolator::set_app_id(int source_app, int target_app)
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

    void Interpolator::compute_interp_coeff(InterpolationMethod /*method = Automatic*/)
    {
        constexpr double eps = 1E-20;
        constexpr double max_d2 = std::numeric_limits<double>::max();

        if (source_bounds_.empty() || target_bounds_.empty())return;

        // TODO: Interpolation method:
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

        node_coeffs_.resize(nn);
        face_coeffs_.resize(nf);

        // init coeff
        for (auto& c : node_coeffs_) { c.src_bd_id = invalid_id; c.dist_sq = max_d2; }
        for (auto& c : face_coeffs_) { c.src_bd_id = invalid_id; c.dist_sq = max_d2; }

        // loop each source bound to accelerate computing.
        for (int ib = 0; ib < static_cast<int>(source_bounds_.size()); ++ib) {
            auto bd_s = source_bounds_.at(ib);

            auto& kdtree = bd_s->kdtree();

            auto it_n = node_coeffs_.begin();
            auto it_f = face_coeffs_.begin();
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
                    if (n == 1 && c.src_bd_id == invalid_id || d2 < c.dist_sq) {
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
                    if (n == 1 && it_f->src_bd_id == invalid_id || d2 < c.dist_sq) {
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

        // do interpolate
        for (int ib = 0; ib < static_cast<int>(source_bounds_.size()); ++ib) {
            auto bd_s = source_bounds_.at(ib);

            auto it_n = node_coeffs_.begin();
            auto it_f = face_coeffs_.begin();
            
            // source does not contain any face: use Local Spline algorithm.
            if (bd_s->face_num() == 0) {
                // loop each node and face-center of target bound
                for (auto bd_t : target_bounds_) {
                    for (auto& p : bd_t->node_coords()) {
                        auto& c = *it_n;

                        // skipping for coincident point or source bound is not this source.
                        if (c.dist_sq <= eps || c.src_bd_id != ib) { ++it_n; continue; }

                        bd_s->compute_local_xps_interp_coeff(
                            p,
                            max_donor,
                            std::span<int_l >(c.donor_nodes,   max_donor),
                            std::span<double>(c.donor_weights, max_donor),
                            c.ndonor
                        );

                        ++it_n;
                    } // next node

                    for (auto& p : bd_t->face_centroids()) {
                        auto& c = *it_f;

                        // skipping for coincident point or source bound is not this source.
                        if (c.dist_sq <= eps || c.src_bd_id != ib) { ++it_f; continue; }

                        bd_s->compute_local_xps_interp_coeff(
                            p,
                            max_donor,
                            std::span<int_l >(c.donor_nodes,   max_donor),
                            std::span<double>(c.donor_weights, max_donor),
                            c.ndonor
                        );

                        ++it_f;
                    } // next face-centroid
                } // next target boundary
            }
            // source contains faces: use projection method
            else {
                // loop each node and face-center of target bound
                int_l ids[8];
                double w[8];
                for (auto bd_t : target_bounds_) {
                    for (auto& p : bd_t->node_coords()) {
                        auto& c = *it_n;

                        // skipping for coincident point or source bound is not this source.
                        if (c.dist_sq <= eps || c.src_bd_id != ib) { ++it_n; continue; }

                        bd_s->compute_project_interp_coeff(
                            p,
                            ids,
                            w,
                            c.ndonor
                        );
                        for (int j = 0; j < c.ndonor; ++j) {
                            c.donor_nodes  [j] = ids[j];
                            c.donor_weights[j] = w  [j];
                        }

                        ++it_n;
                    }// next node

                    for (auto& p : bd_t->face_centroids()) {
                        auto& c = *it_f;

                        // skipping for coincident point or source bound is not this source.
                        if (c.dist_sq <= eps || c.src_bd_id != ib) { ++it_f; continue; }

                        bd_s->compute_project_interp_coeff(
                            p,
                            ids,
                            w,
                            c.ndonor
                        );
                        for (int j = 0; j < c.ndonor; ++j) {
                            c.donor_nodes[j] = ids[j];
                            c.donor_weights[j] = w[j];
                        }

                        ++it_f;
                    } // next face centroid
                } // next target boundary
            }
        }

        computed_ = true;
    }

    void Interpolator::interp_dofs_s2t(std::span<Field* const> sources, std::span<Field*> targets)
    {
        if (!computed_)compute_interp_coeff();
        if (target_bounds_.empty())return;

        // check
        int nerr = 0;
        const int ncomp = targets[0]->info->ncomp;
        for (size_t ib = 0; ib < target_bounds_.size(); ++ib) {
            if (targets[ib]->info->iotype != IncomingDofs) {
                info("***ERROR*** target field \"%s\" is not incoming DOFs!\n", targets[ib]->info->name.c_str());
                ++nerr;
            }
            if (targets[ib]->info->ncomp != ncomp) {
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
            if (sources[ib]->info->ncomp != ncomp) {
                info("***ERROR*** components number of source field \"%s\" not agree!\n", sources[ib]->info->name.c_str());
                ++nerr;
            }
        }
        if (nerr)return;

        int_l nid = 0, fid = 0;
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);
            auto& ft = *(targets[ib]);
            
            // interpolate nodal dofs
            if (ft.info->location == NodeCentered) {
                for (int_l i = 0; i < bt->node_num(); ++i) {
                    auto& coeff = node_coeffs_.at(nid + i);
                    auto& fs = *(sources[coeff.src_bd_id]);
                    auto vals_t = ft.data[i];
                    std::fill(vals_t.begin(), vals_t.end(), 0);
                    for (int j = 0; j < coeff.ndonor; ++j) {
                        auto vals_s = fs.data[coeff.donor_nodes[j]];
                        for (int k = 0; k < ncomp; ++k)
                            vals_t[k] += coeff.donor_weights[k] * vals_s[k];
                    }
                }
            }
            // interpolate face dofs
            else {
                for (int_l i = 0; i < bt->face_num(); ++i) {
                    auto& coeff = face_coeffs_.at(nid + i);
                    auto& fs = *(sources[coeff.src_bd_id]);
                    auto vals_t = ft.data[i];
                    std::fill(vals_t.begin(), vals_t.end(), 0);
                    for (int j = 0; j < coeff.ndonor; ++j) {
                        auto vals_s = fs.data[coeff.donor_nodes[j]];
                        for (int k = 0; k < ncomp; ++k)
                            vals_t[k] += coeff.donor_weights[k] * vals_s[k];
                    }
                }
            }

            nid += bt->node_num();
            fid += bt->face_num();
        }
    }

    void Interpolator::interp_loads_t2s(std::span<Field* const> targets, std::span<Field*> sources)
    {
        if (!computed_)compute_interp_coeff();
        if (source_bounds_.empty())return;

        // check
        int nerr = 0;
        const int ncomp = targets[0]->info->ncomp;
        for (size_t ib = 0; ib < target_bounds_.size(); ++ib) {
            if (targets[ib]->info->iotype != OutgoingLoads) {
                info("***ERROR*** target field \"%s\" is not outgoing LOADs!\n", targets[ib]->info->name.c_str());
                ++nerr;
            }
            if (targets[ib]->info->ncomp != ncomp) {
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
            if (sources[ib]->info->ncomp != ncomp) {
                info("***ERROR*** components number of source field \"%s\" not agree!\n", sources[ib]->info->name.c_str());
                ++nerr;
            }
        }
        if (nerr)return;

        // zero
        for (size_t ib = 0; ib < source_bounds_.size(); ++ib) {
            auto& fs = *(sources[ib]);
            fs.data.fill(0);
        }

        int_l nid = 0, fid = 0;
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);
            auto& ft = *(targets[ib]);

            // interpolate nodal loads
            if (ft.info->location == NodeCentered) {
                for (int_l i = 0; i < bt->node_num(); ++i) {
                    auto& coeff = node_coeffs_.at(nid + i);
                    auto& fs = *(sources[coeff.src_bd_id]);
                    auto vals_t = ft.data[i];
                    for (int j = 0; j < coeff.ndonor; ++j) {
                        auto vals_s = fs.data[coeff.donor_nodes[j]];
                        for (int k = 0; k < ncomp; ++k)
                            vals_s[k] += coeff.donor_weights[k] * vals_t[k];
                    }
                }
            }
            // interpolate face dofs
            else {
                for (int_l i = 0; i < bt->face_num(); ++i) {
                    auto& coeff = face_coeffs_.at(nid + i);
                    auto& fs = *(sources[coeff.src_bd_id]);
                    auto vals_t = ft.data[i];
                    for (int j = 0; j < coeff.ndonor; ++j) {
                        auto vals_s = fs.data[coeff.donor_nodes[j]];
                        for (int k = 0; k < ncomp; ++k)
                            vals_s[k] += coeff.donor_weights[k] * vals_t[k];
                    }
                }
            }

            nid += bt->node_num();
            fid += bt->face_num();
        }
    }
}
