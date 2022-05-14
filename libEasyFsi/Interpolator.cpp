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

    void Interpolator::compute_interp_coeff()
    {
        constexpr double eps = 1E-20;

        if (source_bounds_.empty() || target_bounds_.empty())return;

        // compute global XPS matrix if source boundary does not contain face info.
        //for (Boundary* bd : sources_) {
        //    if (bd->face_num() == 0) {
        //        bd->compute_global_xps_matrix();
        //    }
        //    // use projection method
        //    else {
        //
        //    }
        //}

        int_l nn = 0, nf = 0;
        for (auto bs : source_bounds_) {
            nn += bs->node_num();
            nf += bs->face_num();
        }

        node_coeffs_.resize(nn);
        face_coeffs_.resize(nf);

        int_l nid = 0, fid = 0;
        for (Boundary* bt : target_bounds_) {
            // compute coefficients for each node
            for (int_l i = 0; i < bt->node_num(); ++i, ++ nid) {
                auto& q = bt->node_coords().at(i);
                auto& coeff = node_coeffs_.at(nid);

                // find nearest source boundary
                double d2min = std::numeric_limits<double>::max();
                int isrc = 0, inode = 0;
                for (int ib = 0; ib < static_cast<int>(source_bounds_.size());++ib) {
                    int_l idx = 0;
                    double d2 = 0;
                    source_bounds_.at(ib)->kdtree().search(q.data(), 1, &idx, &d2, 0);
                    if (d2 < d2min) {
                        isrc = ib;
                        inode = idx;
                        d2min = d2;
                    }
                    if (d2 <= eps)break;
                }

                if (d2min <= eps) {
                    coeff.src_bd_id = isrc;
                    coeff.ndonor = 1;
                    coeff.donor_nodes[0] = inode;
                    coeff.donor_weights[0] = 1;
                }
                else {
                    auto bs = source_bounds_.at(isrc);
                    // use projection method
                    if (bs->is_high_order()) {
                        int_l ids[8];
                        double w[8];
                        bs->compute_project_interp_coeff(
                            q,
                            ids,
                            w,
                            coeff.ndonor
                        );
                        for (int j = 0; j < coeff.ndonor; ++j) {
                            coeff.donor_nodes  [j] = ids[j];
                            coeff.donor_weights[j] = w[j];
                        }
                    }
                    // use local spline method
                    else {
                        bs->compute_local_xps_interp_coeff(
                            q,
                            max_donor,
                            std::span<int_l>(coeff.donor_nodes, max_donor),
                            std::span<double>(coeff.donor_weights, max_donor),
                            coeff.ndonor
                        );
                    }
                }
            }

            // compute coefficients for each face
            for (int_l i = 0; i < bt->face_num(); ++i, ++fid) {
                auto& q = bt->face_centroids().at(i);
                auto& coeff = face_coeffs_.at(fid);

                // find nearest source boundary
                double d2min = std::numeric_limits<double>::max();
                int isrc = 0, inode = 0;
                for (int ib = 0; ib < static_cast<int>(source_bounds_.size()); ++ib) {
                    int_l idx = 0;
                    double d2 = 0;
                    source_bounds_.at(ib)->kdtree().search(q.data(), 1, &idx, &d2, 0);
                    if (d2 < d2min) {
                        isrc = ib;
                        inode = idx;
                        d2min = d2;
                    }
                    if (d2 <= eps)break;
                }

                if (d2min <= eps) {
                    coeff.src_bd_id = isrc;
                    coeff.ndonor = 1;
                    coeff.donor_nodes[0] = inode;
                    coeff.donor_weights[0] = 1;
                }
                else {
                    auto bs = source_bounds_.at(isrc);
                    // use projection method
                    int_l ids[8];
                    double w[8];
                    bs->compute_project_interp_coeff(
                        q,
                        ids,
                        w,
                        coeff.ndonor
                    );
                    for (int j = 0; j < coeff.ndonor; ++j) {
                        coeff.donor_nodes[j] = ids[j];
                        coeff.donor_weights[j] = w[j];
                    }
                }
            }
        }

        computed_ = true;
    }

    void Interpolator::interp_target_dofs(std::span<Field* const> sources, std::span<Field*> targets)
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

    void Interpolator::interp_source_loads(std::span<Field* const> targets, std::span<Field*> sources)
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
