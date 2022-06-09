#include <fstream>

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
        node_info_.clear();
        face_info_.clear();
        //donors_.clear();
        //weights_.clear();
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

        // do interpolate
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
        if (!computed_)error("coefficients are not computed!");
        if (target_bounds_.empty())return;

        // check
        int nerr = 0;
        const int ncomp = targets[0]->info->ncomp;
        for (size_t ib = 0; ib < target_bounds_.size(); ++ib) {
            if (targets[ib] == nullptr)continue;

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

        auto nit = node_info_.begin();
        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);

            if (targets[ib] == nullptr) {
                nit += bt->node_num();
                fit += bt->face_num();
                continue;
            }

            auto& ft = *(targets[ib]);
            
            // interpolate nodal dofs
            if (ft.info->location == NodeCentered) {
                for (int_l i = 0; i < bt->node_num(); ++i, ++nit) {
                    auto& coeff = *nit;
                    auto& fs = *(sources[coeff.src_bd_id]); // source field
                    auto vals_t = ft.data[i];
                    std::fill(vals_t.begin(), vals_t.end(), 0);
                    for (int j = 0; j < coeff.ndonor; ++j) {
                        auto donor    = coeff.donor_nodes[j];
                        auto weight   = coeff.donor_weights[j];
                        auto vals_s   = fs.data[donor];
                        for (int k = 0; k < ncomp; ++k)
                            vals_t[k] += weight * vals_s[k];
                    }
                }
                fit += bt->face_num();
            }
            // interpolate face dofs
            else {
                for (int_l i = 0; i < bt->face_num(); ++i, ++fit) {
                    auto& coeff = *fit;
                    auto& fs = *(sources[coeff.src_bd_id]); // source field
                    auto vals_t = ft.data[i];
                    std::fill(vals_t.begin(), vals_t.end(), 0);
                    for (int j = 0; j < coeff.ndonor; ++j) {
                        auto donor = coeff.donor_nodes[j];
                        auto weight = coeff.donor_weights[j];
                        auto vals_s = fs.data[donor];
                        for (int k = 0; k < ncomp; ++k)
                            vals_t[k] += weight * vals_s[k];
                    }
                }
                nit += bt->node_num();
            }
        }
    }

    void Interpolator::interp_load_t2s(std::span<Field* const> targets, std::span<Field*> sources, bool fill_src_zeros_first/* = true*/)const
    {
        if (!computed_)error("coefficients are not computed!");
        if (source_bounds_.empty())return;

        // check
        int nerr = 0;
        const int ncomp = targets[0]->info->ncomp;
        for (size_t ib = 0; ib < target_bounds_.size(); ++ib) {
            if (targets[ib] == nullptr)continue;

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
        if (fill_src_zeros_first) {
            for (size_t ib = 0; ib < source_bounds_.size(); ++ib) {
                auto& fs = *(sources[ib]);
                fs.data.fill(0);
            }
        }

        auto nit = node_info_.begin();
        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);
            if (targets[ib] == nullptr) {
                nit += bt->node_num();
                fit += bt->face_num();
                continue;
            }

            auto& ft = *(targets[ib]);

            // interpolate nodal loads
            if (ft.info->location == NodeCentered) {
                for (int_l i = 0; i < bt->node_num(); ++i, ++nit) {
                    auto& coeff = *nit;
                    auto& fs = *(sources[coeff.src_bd_id]); // source field
                    auto vals_t = ft.data[i];
                    for (int j = 0; j < coeff.ndonor; ++j) {
                        auto donor = coeff.donor_nodes[j];
                        auto weight = coeff.donor_weights[j];
                        auto vals_s = fs.data[donor];
                        for (int k = 0; k < ncomp; ++k)
                            vals_s[k] += weight * vals_t[k];
                    }
                }
                fit += bt->face_num();
            }
            // interpolate face dofs
            else {
                for (int_l i = 0; i < bt->face_num(); ++i, ++fit) {
                    auto& coeff = *fit;
                    auto& fs = *(sources[coeff.src_bd_id]); // source field
                    auto vals_t = ft.data[i];
                    for (int j = 0; j < coeff.ndonor; ++j) {
                        auto donor = coeff.donor_nodes[j];
                        auto weight = coeff.donor_weights[j];
                        auto vals_s = fs.data[donor];
                        for (int k = 0; k < ncomp; ++k)
                            vals_s[k] += weight * vals_t[k];
                    }
                }
                nit += bt->node_num();
            }
        }
    }

    void Interpolator::interp_node_dofs_s2t(int nfield, const double** src_node_dofs, double** des_node_dofs)const
    {
        auto nit = node_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);
            if (des_node_dofs[ib] == nullptr) {
                nit += bt->node_num();
                continue;
            }

            auto vals_t = des_node_dofs[ib];

            // interpolate nodal dofs
            for (int_l i = 0; i < bt->node_num(); ++i, ++nit) {
                auto& coeff = *nit;
                auto  fs = src_node_dofs[coeff.src_bd_id]; // source field
                std::fill(vals_t, vals_t + nfield, 0);
                for (int j = 0; j < coeff.ndonor; ++j) {
                    auto donor = coeff.donor_nodes[j];
                    auto weight = coeff.donor_weights[j];
                    auto vals_s = fs + nfield * donor;
                    for (int k = 0; k < nfield; ++k, ++vals_t)
                        vals_t[k] += weight * vals_s[k];
                }
                vals_t += nfield;
            }
        }
    }
    void Interpolator::interp_face_dofs_s2t(int nfield, const double** src_node_dofs, double** des_face_dofs)const
    {
        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);
            if (des_face_dofs[ib] == nullptr) {
                fit += bt->face_num();
                continue;
            }

            auto vals_t = des_face_dofs[ib];

            // interpolate face dofs
            for (int_l i = 0; i < bt->face_num(); ++i, ++fit) {
                auto& coeff = *fit;
                auto  fs = src_node_dofs[coeff.src_bd_id]; // source field
                std::fill(vals_t, vals_t + nfield, 0);
                for (int j = 0; j < coeff.ndonor; ++j) {
                    auto donor = coeff.donor_nodes[j];
                    auto weight = coeff.donor_weights[j];
                    auto vals_s = fs + nfield * donor;
                    for (int k = 0; k < nfield; ++k, ++vals_t)
                        vals_t[k] += weight * vals_s[k];
                }
                vals_t += nfield;
            }
        }
    }
    void Interpolator::interp_node_load_t2s(int nfield, double** src_node_load, const double** des_node_load, bool fill_src_zeros_first/* = true*/)const
    {
        if (fill_src_zeros_first) {
            for (size_t ib = 0; ib < source_bounds_.size(); ++ib) {
                std::memset(src_node_load[ib], 0, sizeof(double) * nfield * source_bounds_.at(ib)->node_num());
            }
        }

        auto nit = node_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);
            if (des_node_load[ib] == nullptr) {
                nit += bt->node_num();
                continue;
            }

            auto vals_t = des_node_load[ib];

            // interpolate nodal loads
            for (int_l i = 0; i < bt->node_num(); ++i, ++nit) {
                auto& coeff = *nit;
                auto fs = src_node_load[coeff.src_bd_id]; // source field
                for (int j = 0; j < coeff.ndonor; ++j) {
                    auto donor = coeff.donor_nodes[j];
                    auto weight = coeff.donor_weights[j];
                    auto vals_s = fs + donor * nfield;
                    for (int k = 0; k < nfield; ++k)
                        vals_s[k] += weight * vals_t[k];
                }
                vals_t += nfield;
            }
        }
    }
    void Interpolator::interp_face_load_t2s(int nfield, double** src_node_load, const double** des_face_load, bool fill_src_zeros_first/* = true*/)const
    {
        if (fill_src_zeros_first) {
            for (size_t ib = 0; ib < source_bounds_.size(); ++ib) {
                std::memset(src_node_load[ib], 0, sizeof(double) * nfield * source_bounds_.at(ib)->node_num());
            }
        }

        auto fit = face_info_.begin();
        for (int ib = 0; ib < target_bounds_.size(); ++ib) {
            auto& bt = target_bounds_.at(ib);
            if (des_face_load[ib] == nullptr) {
                fit += bt->face_num();
                continue;
            }

            auto vals_t = des_face_load[ib];

            // interpolate nodal loads
            for (int_l i = 0; i < bt->face_num(); ++i, ++fit) {
                auto& coeff = *fit;
                auto fs = src_node_load[coeff.src_bd_id]; // source field
                for (int j = 0; j < coeff.ndonor; ++j) {
                    auto donor = coeff.donor_nodes[j];
                    auto weight = coeff.donor_weights[j];
                    auto vals_s = fs + donor * nfield;
                    for (int k = 0; k < nfield; ++k)
                        vals_s[k] += weight * vals_t[k];
                }
                vals_t += nfield;
            }
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

    void Interpolator::read_coefficients(const char* file)
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

    void Interpolator::write_tecplot(const char* file, int nfield, const char** field_names, const double** src_node_fields, const double** des_node_fields)const
    {
        std::ofstream fout(file);
        if (!fout) {
            error("%s(), failed open file: %s\n", __func__, file);
            return;
        }

        info("\nwriting Tecplot file: %s\n", file);

        int izone = 0;
        int sp_zone = 0, fp_zone = 0;

        //--- file header

        fout
            << "TITLE = \"Interpolator Data for FSI\"\n"
            << "FILETYPE = FULL\n"
            << "VARIABLES = \"X\",\"Y\",\"Z\"";
        for (int i = 0; i < nfield; ++i) {
            fout << ",\"" << field_names[i] << "\"";
        }
        fout << '\n';

        //--- solid zones

        int iz = 0;
        for (auto bd : source_bounds_) {
            int nns = bd->node_num();

            // points zone
            {
                auto pnt_s = bd->node_coords().data();
                fout << "ZONE T=\""<< bd->name() << ".Nodes\" ZONETYPE=ORDERED I=" << nns << " DATAPACKING = POINT\n";
                for (int i = 0; i < nns; ++i) {
                    fout
                        << std::scientific << pnt_s[i].x << ' '
                        << std::scientific << pnt_s[i].y << ' '
                        << std::scientific << pnt_s[i].z;
                    for (int j = 0; j < nfield; ++j) {
                        fout << ' ' << std::scientific << src_node_fields[iz][i * nfield + j];
                    }
                    fout << '\n';
                }
                ++izone;
                sp_zone = izone;
            }

            // elements zone
            if      (bd->topo() == ZT_CURVE) {
                fout
                    << "ZONE T=\"" << bd->name() << ".Elements\""
                    << " ZONETYPE=FELINESEG"
                    << " NODES=" << nns
                    << " ELEMENTS=" << bd->face_num()
                    << " VARSHARELIST=([1-" << 3 + nfield << "]=" << sp_zone << ")\n";

                auto& fnodes = bd->face_nodes();
                for (int i = 0; i < bd->face_num(); ++i) {
                    fout
                        << fnodes(i, 0) + 1 << ' '
                        << fnodes(i, 1) + 1 << '\n';
                }
                ++izone;
            }
            else if (bd->topo() == ZT_SURFACE) {

                auto& fts = bd->face_types();
                auto& fnodes = bd->face_nodes();

                // all are triangles
                if (std::all_of(fts.begin(), fts.end(), [](int ft) {return ft == FT_TRI3 || ft == FT_TRI6; })) {
                    fout
                        << "ZONE T=\"" << bd->name() << ".Elements\""
                        << " ZONETYPE=FETRIANGLE"
                        << " NODES=" << nns
                        << " ELEMENTS=" << bd->face_num()
                        << " VARSHARELIST=([1-" << 3 + nfield << "]=1)\n";

                    for (int i = 0; i < bd->face_num(); ++i) {
                        fout
                            << fnodes(i, 0) + 1 << ' '
                            << fnodes(i, 1) + 1 << ' '
                            << fnodes(i, 2) + 1 << '\n';
                    }
                    ++izone;
                }
                // mixed elements, save as quad
                else if (!bd->contains_polygon()) {
                    fout
                        << "ZONE T=\"" << bd->name() << ".Elements\""
                        << " ZONETYPE=FEQUADRILATERAL"
                        << " NODES=" << nns
                        << " ELEMENTS=" << bd->face_num()
                        << " VARSHARELIST=([1-" << 3 + nfield << "]=" << sp_zone << ")\n";

                    for (int i = 0; i < bd->face_num(); ++i) {
                        if (fts[i] == FT_TRI3 || fts[i] == FT_TRI6) {
                            fout
                                << fnodes(i, 0) + 1 << ' '
                                << fnodes(i, 1) + 1 << ' '
                                << fnodes(i, 2) + 1 << ' '
                                << fnodes(i, 0) + 1 << '\n';
                        }
                        else {
                            fout
                                << fnodes(i, 0) + 1 << ' '
                                << fnodes(i, 1) + 1 << ' '
                                << fnodes(i, 2) + 1 << ' '
                                << fnodes(i, 3) + 1 << '\n';
                        }
                    }
                    ++izone;
                }
                // polygon zone
                else {
                    fout
                        << "ZONE T=\"Solid.Elements\""
                        << " ZONETYPE=FEPOLYGON"
                        << " NODES=" << nns
                        << " FACES=" << bd->edges_for_surface().size()
                        << " ELEMENTS=" << bd->face_num()
                        << " VARSHARELIST=([1-" << 3 + nfield << "]=1)\n";

                    // face nodes
                    for (auto& e : bd->edges_for_surface()) {
                        fout
                            << e.n0 + 1 << ' '
                            << e.n1 + 1 << '\n';
                    }

                    // left elements
                    int j = 1;
                    for (auto& e : bd->edges_for_surface()) {
                        fout << e.f0 + 1 << ((j % 8 != 0) ? ' ' : '\n');
                        ++j;
                    }
                    fout << '\n';

                    // right elements (negative indicates boundary connections)
                    j = 1;
                    for (auto& e : bd->edges_for_surface()) {
                        fout << (e.f1 != invalid_id ? e.f1 + 1 : -1) << ((j % 8 != 0) ? ' ' : '\n');
                        ++j;
                    }
                    fout << '\n';

                    ++izone;
                }
            }
            else if (bd->topo() == ZT_VOLUME) {
                // TBD
                //++izone;
                error("volume zone is not implemented");
            }

            ++iz;
        }
        
        //--- fluid zone

        iz = 0;
        for (auto bd : target_bounds_) {
            // points zone
            int nnf = bd->node_num();
            int nff = bd->face_num();
            {
                auto pnt_f = bd->node_coords().data();
                fout << "ZONE T=\"" << bd->name() << ".Nodes\" ZONETYPE=ORDERED I=" << nnf << " DATAPACKING = POINT\n";
                for (int i = 0; i < nnf; ++i) {
                    fout
                        << std::scientific << pnt_f[i].x << ' '
                        << std::scientific << pnt_f[i].y << ' '
                        << std::scientific << pnt_f[i].z;
                    for (int j = 0; j < nfield; ++j) {
                        fout << ' ' << std::scientific << des_node_fields[iz][i * nfield + j];
                    }
                    fout << '\n';
                }
                ++izone;
                fp_zone = izone;
            }

            //--- elements zone

            if      (bd->topo() == ZT_CURVE) {
                fout
                    << "ZONE T=\"" << bd->name() << ".Elements\""
                    << " ZONETYPE=FELINESEG"
                    << " NODES=" << nnf
                    << " ELEMENTS=" << nff
                    << " VARSHARELIST=([1-" << 3 + nfield << "] = " << fp_zone << ")\n";

                auto& fnodes = bd->face_nodes();
                for (int i = 0; i < nff; ++i) {
                    fout
                        << fnodes(i, 0) + 1 << ' '
                        << fnodes(i, 1) + 1 << '\n';
                }
                ++izone;
            }
            else if (bd->topo() == ZT_SURFACE) {

                auto& fts    = bd->face_types();
                auto& fnodes = bd->face_nodes();

                // all faces are triangle
                if (std::all_of(fts.begin(), fts.end(), [](int ft) {return ft == FT_TRI3 || ft == FT_TRI6; })) {
                    fout
                        << "ZONE T=\"" << bd->name() << ".Elements\""
                        << " ZONETYPE=FETRIANGLE"
                        << " NODES=" << nnf
                        << " ELEMENTS=" << bd->face_num()
                        << " VARSHARELIST=([1-" << 3 + nfield << "] = " << fp_zone << ")\n";

                    for (int i = 0; i < nff; ++i) {
                        fout
                            << fnodes(i, 0) + 1 << ' '
                            << fnodes(i, 1) + 1 << ' '
                            << fnodes(i, 2) + 1 << '\n';
                    }
                    ++izone;
                }
                // not polygon
                else if (!bd->contains_polygon()) {
                    fout
                        << "ZONE T=\"" << bd->name() << ".Elements\""
                        << " ZONETYPE=FEQUADRILATERAL"
                        << " NODES=" << nnf
                        << " ELEMENTS=" << nff
                        << " VARSHARELIST=([1-" << 3 + nfield << "] = " << fp_zone << ")\n";

                    for (int i = 0; i < nff; ++i) {
                        if (fts[i] == FT_TRI3 || fts[i] == FT_TRI6) {
                            fout
                                << fnodes(i, 0) + 1 << ' '
                                << fnodes(i, 1) + 1 << ' '
                                << fnodes(i, 2) + 1 << ' '
                                << fnodes(i, 0) + 1 << '\n';
                        }
                        else {
                            fout
                                << fnodes(i, 0) + 1 << ' '
                                << fnodes(i, 1) + 1 << ' '
                                << fnodes(i, 2) + 1 << ' '
                                << fnodes(i, 3) + 1 << '\n';
                        }
                    }
                    ++izone;
                }
                // polygon zone
                else {
                    fout
                        << "ZONE T=\"" << bd->name() << ".Elements\""
                        << " ZONETYPE=FEQUADRILATERAL"
                        << " NODES=" << bd->node_num()
                        << " FACES=" << bd->edges_for_surface().size()
                        << " ELEMENTS=" << bd->face_num()
                        << " VARSHARELIST=([1-" << 3 + nfield << "] = " << fp_zone << ")\n";

                    // face nodes
                    for (auto& e : bd->edges_for_surface()) {
                        fout
                            << e.n0 + 1 << ' '
                            << e.n1 + 1 << '\n';
                    }

                    // left elements
                    int j = 1;
                    for (auto& e : bd->edges_for_surface()) {
                        fout << e.f0 + 1 << ((j % 8 != 0) ? ' ' : '\n');
                        ++j;
                    }
                    fout << '\n';

                    // right elements (negative indicates boundary connections)
                    j = 1;
                    for (auto& e : bd->edges_for_surface()) {
                        fout << (e.f1 != invalid_id ? e.f1 + 1 : -1) << ((j % 8 != 0) ? ' ' : '\n');
                        ++j;
                    }
                    fout << '\n';

                    ++izone;
                }
            }
            else if (bd->topo() == ZT_VOLUME) {
                // TBD
                //++izone;
            }

            ++iz;
        }
        

        fout.close();

        info("!!!OK!!!\n");
    }
}
