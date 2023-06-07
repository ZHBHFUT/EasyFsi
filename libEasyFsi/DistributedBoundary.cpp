#include "Logger.hpp"
#include "DistributedBoundary.hpp"

namespace EasyLib {

    void DistributedBoundary::clear()
    {
        full_bound_.clear();

        part_nodes_ia_g_.clear(); part_nodes_ia_g_.push_back(0);
        part_nodes_ja_g_.clear();
        part_faces_ia_g_.clear(); part_faces_ia_g_.push_back(0);

        local_bound_ = nullptr;
        comm_ = nullptr;
        root_ = 0;
    }

    void DistributedBoundary::assemble(Boundary& local_bound, Communicator& comm, int root_rank)
    {
        clear();

        comm_ = &comm;
        root_ = root_rank;
        local_bound_ = &local_bound;

        if (comm_->size() <= 1) {
            full_bound_ = local_bound;
            part_nodes_ia_g_.push_back(local_bound.node_num());
            part_faces_ia_g_.push_back(local_bound.face_num());
            part_nodes_ja_g_.resize(local_bound.node_num());
            for (int_l i = 0; i < local_bound.node_num(); ++i)
                part_nodes_ja_g_[i] = i;
            //info("\nassemble boundary: do nothing for serial running.\n");
            return;
        }

        info("\nassemble boundary\n");

        constexpr int TAG_BD_SIZE = 0;
        constexpr int TAG_BD_DATA = 100;

        int_l buf[3];
        if (comm_->rank() == root_) {

            // copy fields definition
            //for (auto& f : local_bound_->get_fields()) {
            //    full_bound_.register_field(
            //        f.name.c_str(),
            //        f.ncomp,
            //        f.location,
            //        f.iotype,
            //        f.units.c_str()
            //    );
            //}

            // name
            full_bound_.set_name(local_bound.name().c_str());
            full_bound_.set_user_id(local_bound.user_id());

            part_nodes_ia_g_.resize(comm_->size() + 1, 0);
            part_faces_ia_g_.resize(comm_->size() + 1, 0);

            int_l nface_g = 0, nfnodes_g = 0, max_nnode_g = 0;
            for (int part = 0; part < comm_->size(); ++part) {
                if (part != root_rank) {
                    comm_->recv(buf, 3, part, TAG_BD_SIZE);
                }
                else {
                    buf[0] = local_bound.node_num();
                    buf[1] = local_bound.face_num();
                    buf[2] = local_bound.face_nodes().ndata();
                }
                max_nnode_g += buf[0];
                nface_g     += buf[1];
                nfnodes_g   += buf[2];
                part_nodes_ia_g_[part + 1] += buf[0];
                part_faces_ia_g_[part + 1] += buf[1];
            }

            // allocate
            full_bound_.reserve(max_nnode_g, nface_g, nfnodes_g);
            part_nodes_ja_g_.resize(max_nnode_g);

            Boundary bnd;
            ivec fnodes;
            for (int part = 0; part < comm_->size(); ++part) {
                const Boundary* pb;
                if (part != root_rank) {
                    comm_->recv(bnd, part, TAG_BD_DATA);
                    pb = &bnd;
                }
                else
                    pb = &local_bound;

                // add nodes to set
                const auto* coords = reinterpret_cast<const vec3*>(pb->node_coords().data());
                for (int_l i = 0; i < pb->node_num(); ++i) {
                    part_nodes_ja_g_[part_nodes_ia_g_[part] + i] = full_bound_.add_node(coords[i], pb->nodes().l2g(i));
                }

                // add faces
                const auto* fcents = reinterpret_cast<const vec3*>(pb->face_centroids().data());
                for (int_l i = 0; i < pb->face_num(); ++i) {
                    auto nodes = pb->face_nodes()[i];
                    fnodes.resize(nodes.size());
                    for (size_t j = 0; j < nodes.size(); ++j) {
                        auto g = pb->nodes().l2g(nodes[j]);
                        fnodes[j] = full_bound_.nodes().g2l(g);
                    }
                    full_bound_.add_face((FaceTopo)pb->face_types()[i], static_cast<int_l>(nodes.size()), fnodes.data(), fcents[i]);
                }

                //TODO: check fields
            }

            info(
                "    Total Node = %d\n"
                "    Total Face = %d\n",
                full_bound_.node_num(),
                full_bound_.face_num()
            );

            // allocate fields
            //for (auto& f : full_bound_.get_fields()) {
            //    if (f.location == NodeCentered)
            //        f.data.resize(full_bound_.node_num(), f.ncomp);
            //    else
            //        f.data.resize(full_bound_.face_num(), f.ncomp);
            //}

            // compute metrics
            full_bound_.compute_metics();
        }
        else {
            // send size info
            buf[0] = local_bound.node_num();
            buf[1] = local_bound.face_num();
            buf[2] = local_bound.face_nodes().ndata();
            comm_->send(buf, 3, root_, TAG_BD_SIZE);
            comm_->send(local_bound, root_, TAG_BD_DATA);
        }
    }

    void DistributedBoundary::gather_node_fields(int nfields, const double* local_fields, double* global_fields)
    {
        //--- serial running

        const int_l count = local_bound_ ? local_bound_->node_num() : 0;
        if (!comm_ || comm_->size() <= 1) {
            if (local_bound_)
                std::memcpy(global_fields, local_fields, sizeof(double) * nfields * count);
            return;
        }

        //--- parallel running

        if (comm_->rank() == root_) {
            buffer_g_.resize(nfields * full_bound_.node_num());
            auto& ia = part_nodes_ia_g_;
            auto& ja = part_nodes_ja_g_;

            for (int ip = 0; ip < comm_->size(); ++ip) {
                auto offset = ia[ip];
                auto n = ia[ip + 1] - offset; // node number
                if (n <= 0)continue;

                // receive data
                if (ip != root_)
                    comm_->recv(buffer_g_.data(), nfields * n, ip, 0);

                auto buf = ip != root_ ? buffer_g_.data() : local_fields;

                // assemble global data
                for (int i = 0; i < n; ++i) {
                    auto id = ja[offset + i];
                    for (int j = 0; j < nfields; ++j) {
                        global_fields[nfields * id + j] = buf[nfields * i + j];
                    }
                }
            }
        }
        else if (count > 0) {
            comm_->send(local_fields, nfields * count, root_, 0);
        }
    }
    void DistributedBoundary::gather_face_fields(int nfields, const double* local_fields, double* global_fields)
    {
        //--- serial running

        const int_l count = local_bound_ ? local_bound_->face_num() : 0;
        if (!comm_ || comm_->size() <= 1) {
            std::memcpy(global_fields, local_fields, sizeof(double) * nfields * count);
            return;
        }

        //--- parallel running

        if (comm_->rank() == root_) {
            auto& ia = part_faces_ia_g_;
            for (int ip = 0; ip < comm_->size(); ++ip) {
                auto offset = ia[ip];
                auto n = ia[ip + 1] - offset; // face number
                if (n <= 0)continue;

                // receive data
                if (ip != root_) {
                    comm_->recv(global_fields + offset * nfields, nfields * n, ip, 0);
                }
                // directly copy
                else {
                    memcpy(global_fields + offset * nfields, local_fields, sizeof(double) * nfields * n);
                }
            }
        }
        else if (count > 0) {
            comm_->send(local_fields, nfields * count, root_, 0);
        }
    }
    void DistributedBoundary::accumulate_node_fields(int nfields, const double* local_fields, double* global_fields)
    {
        //--- serial running

        const int_l count = local_bound_ ? local_bound_->node_num() : 0;
        if (!comm_ || comm_->size() <= 1) {
            std::memcpy(global_fields, local_fields, sizeof(double) * nfields * count);
            return;
        }

        //--- parallel running

        if (comm_->rank() == root_) {
            buffer_g_.resize(nfields * full_bound_.node_num());
            auto& ia = part_nodes_ia_g_;
            auto& ja = part_nodes_ja_g_;

            for (int ip = 0; ip < comm_->size(); ++ip) {
                auto offset = ia[ip];
                auto n = ia[ip + 1] - offset; // node number
                if (n <= 0)continue;

                // receive data
                if (ip != root_)
                    comm_->recv(buffer_g_.data(), nfields * n, ip, 0);

                auto buf = ip != root_ ? buffer_g_.data() : local_fields;

                // assemble global data
                for (int i = 0; i < n; ++i) {
                    auto id = ja[offset + i];
                    for (int j = 0; j < nfields; ++j) {
                        global_fields[nfields * id + j] += buf[nfields * i + j];
                    }
                }
            }
        }
        else if (count > 0) {
            comm_->send(local_fields, nfields * count, root_, 0);
        }
    }
    void DistributedBoundary::accumulate_face_fields(int nfields, const double* local_fields, double* global_fields)
    {
        //--- serial running

        const int_l count = local_bound_ ? local_bound_->face_num() : 0;
        if (!comm_ || comm_->size() <= 1) {
            for (int i = 0; i < count; ++i) {
                for (int j = 0; j < nfields; ++j) {
                    global_fields[nfields * i + j] = local_fields[nfields * i + j];
                }
            }
            return;
        }

        //--- parallel running

        if (comm_->rank() == root_) {
            buffer_g_.resize(nfields * full_bound_.face_num());
            auto& ia = part_faces_ia_g_;
            for (int ip = 0; ip < comm_->size(); ++ip) {
                auto offset = ia[ip];
                auto n = ia[ip + 1] - offset; // face number
                if (n <= 0)continue;

                // receive data
                if (ip != root_)
                    comm_->recv(buffer_g_.data(), nfields * n, ip, 0);

                auto src = ip != root_ ? buffer_g_.data() : local_fields;
                auto des = global_fields + offset * nfields;
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < nfields; ++j) {
                        des[nfields * i + j] += src[nfields * i + j];
                    }
                }
            }
        }
        else if (count > 0) {
            comm_->send(local_fields, nfields * count, root_, 0);
        }
    }
    void DistributedBoundary::scatter_node_fields(int nfields, const double* global_fields, double* local_fields)
    {
        //--- serial running

        const int_l count = local_bound_ ? local_bound_->node_num() : 0;
        if (!comm_ || comm_->size() <= 1) {
            std::memcpy(local_fields, global_fields, sizeof(double) * nfields * count);
            return;
        }

        //--- parallel running

        if (comm_->rank() == root_) {
            buffer_g_.resize(nfields * full_bound_.node_num());
            auto& ia = part_nodes_ia_g_;
            auto& ja = part_nodes_ja_g_;

            for (int ip = 0; ip < comm_->size(); ++ip) {
                auto offset = ia[ip];
                auto n = ia[ip + 1] - offset; // node number
                if (n <= 0)continue;

                auto buf = ip != root_ ? buffer_g_.data() : local_fields;

                // copy to buffer
                for (int i = 0; i < n; ++i) {
                    auto id = ja[offset + i];
                    for (int j = 0; j < nfields; ++j) {
                        buf[nfields * i + j] = global_fields[nfields * id + j];
                    }
                }

                // send data
                if (ip != root_)
                    comm_->send(buffer_g_.data(), nfields * n, ip, 0);
            }
        }
        else if (count > 0) {
            comm_->recv(local_fields, nfields * count, root_, 0);
        }
    }
    void DistributedBoundary::scatter_face_fields(int nfields, const double* global_fields, double* local_fields)
    {
        //--- serial running

        const int_l count = local_bound_ ? local_bound_->face_num() : 0;
        if (!comm_ || comm_->size() <= 1) {
            std::memcpy(local_fields, global_fields, sizeof(double) * nfields * count);
            return;
        }

        //--- parallel running

        if (comm_->rank() == root_) {
            auto& ia = part_faces_ia_g_;
            for (int ip = 0; ip < comm_->size(); ++ip) {
                auto offset = ia[ip];
                auto n = ia[ip + 1] - offset; // face number
                if (n <= 0)continue;

                // send data
                if (ip != root_) {
                    comm_->send(global_fields + offset * nfields, nfields * n, ip, 0);
                }
                // directly copy
                else {
                    memcpy(local_fields, global_fields + offset * nfields, sizeof(double) * nfields * n);
                }
            }
        }
        else if (count > 0) {
            comm_->recv(local_fields, nfields * count, root_, 0);
        }
    }

    //void DistributedBoundary::gather_fields(const char* field_name)
    //{
    //    double* local_fields = nullptr;
    //    double* global_fields = nullptr;
    //
    //    auto& f = local_bound_->get_field(field_name);
    //    local_fields = f.data.data();
    //
    //    if (comm_->rank() == root_) {
    //        auto& gf = full_bound_.get_field(field_name);
    //        global_fields = gf.data.data();
    //    }
    //
    //    if (f.location == NodeCentered) {
    //        if (f.iotype == IncomingDofs || f.iotype == OutgoingDofs)
    //            gather_node_fields(f.ncomp, local_fields, global_fields);
    //        else
    //            accumulate_node_fields(f.ncomp, local_fields, global_fields);
    //    }
    //    else {
    //        gather_face_fields(f.ncomp, local_fields, global_fields);
    //    }
    //}
    //void DistributedBoundary::scatter_fields(const char* field_name)
    //{
    //    double* local_fields = nullptr;
    //    double* global_fields = nullptr;
    //
    //    auto& f = local_bound_->get_field(field_name);
    //    local_fields = f.data.data();
    //
    //    if (comm_->rank() == root_) {
    //        auto& gf = full_bound_.get_field(field_name);
    //        global_fields = gf.data.data();
    //    }
    //
    //    if (f.location == NodeCentered) {
    //        scatter_node_fields(f.ncomp, global_fields, local_fields);
    //    }
    //    else {
    //        scatter_face_fields(f.ncomp, global_fields, local_fields);
    //    }
    //}
}
