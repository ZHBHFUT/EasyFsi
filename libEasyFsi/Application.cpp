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
//! @file       Application.cpp
//!             The implement of Application class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <fstream>

#include "Logger.hpp"
#include "Application.hpp"

namespace EasyLib {

    Application::Application(const char* name)
    {
        data_.app_name = name;
    }

    Application::Application(const char* name, Communicator& intra_comm, int root)
    {
        data_.app_name = name;
        intra_comm_ = &intra_comm;
        intra_root_ = root;
    }

    void Application::clear()
    {
        stop_coupling();

        inter_comm_ = nullptr;
        intra_comm_ = nullptr;
        intra_root_ = 0;

        data_.bounds.clear();
        data_.field_defs.clear();
        data_.iter = 0;
        data_.rank = 0;
        data_.time = 0;

        remote_apps_.clear();
        local_bounds_.clear();
        bounds_.clear();
        is_started_ = false;
        interps_.clear();
        field_interps_.clear();

        getter_ = nullptr;
        setter_ = nullptr;
        this_fields_.clear();
        remote_fields_.clear();
    }

    //void Application::create(const char* name)
    //{
    //    clear();
    //    data_.app_name = name;
    //}
    //
    //void Application::create(const char* name, Communicator& intra_comm, int root)
    //{
    //    clear();
    //    data_.app_name = name;
    //    intra_comm_ = &intra_comm;
    //    intra_root_ = root;
    //}

    Boundary& Application::add_coupled_boundary()
    {
        if (is_started_) {
            if (!intra_comm_ || intra_comm_->rank() == intra_root_)
                warn("unable add boundary during solving!");
            return *(Boundary*)nullptr;
        }
#if __cplusplus >= 201703L
        return local_bounds_.emplace_back(Boundary{});
#else
        local_bounds_.push_back(Boundary{});
        return local_bounds_.back();
#endif
    }

    void Application::set_field_function(get_boundary_field_function getter, set_boundary_field_function setter)
    {
        getter_ = getter;
        setter_ = setter;
    }

    void Application::register_field(const char* field_name, int ncomp, FieldLocation location, FieldIO iotype, const char* units)
    {
        if (is_started_) {
            if (!intra_comm_ || intra_comm_->rank() == intra_root_)
                warn("unable register new field during solving!");
            return;
        }

        if (!field_name || *field_name == '\0') {
            if (!intra_comm_ || intra_comm_->rank() == intra_root_)
                warn("invalid field name!");
            return;
        }
        if (ncomp <= 0) {
            if (!intra_comm_ || intra_comm_->rank() == intra_root_)
                warn("invalid component number of field!");
            return;
        }

        for (auto& f : data_.field_defs) {
            if (f.name == field_name) {
                if (!intra_comm_ || intra_comm_->rank() == intra_root_)
                    warn("field \"%s\" already registered!", field_name);
                f.ncomp = ncomp;
                f.location = location;
                f.iotype = iotype;
                f.units = units ? units : "";
                return;
            }
        }

#if __cplusplus >= 201703L
        auto& f = data_.field_defs.emplace_back(FieldInfo{});
#else
        data_.field_defs.push_back(FieldInfo{});
        auto& f = data_.field_defs.back();
#endif
        f.name = field_name;
        f.ncomp = ncomp;
        f.location = location;
        f.iotype = iotype;
        f.units = units ? units : "";
    }

    bool Application::start_coupling(Communicator& comm)
    {
        inter_comm_ = &comm;

        if (!intra_comm_ || intra_comm_->rank() == intra_root_)
            info("\nstart coupling!\n");

        // assemble boundary if has intra-communicator, and use assembled boundary
        if (intra_comm_ && intra_comm_->size() > 1) {
            bounds_.clear();
            bounds_.resize(local_bounds_.size());
            data_.bounds.resize(local_bounds_.size());
            for (size_t i = 0; i < bounds_.size(); ++i) {
                auto& b = bounds_.at(i);
                b.assemble(local_bounds_.at(i), *intra_comm_, intra_root_);
                data_.bounds.at(i) = &b.full_boundary();
            }
        }
        // use local boundary
        else {
            bounds_.clear();
            data_.bounds.resize(local_bounds_.size());
            for (size_t i = 0; i < local_bounds_.size(); ++i) {
                data_.bounds.at(i) = &local_bounds_.at(i);
            }
        }

        // clear remote apps.
        for (auto& app : remote_apps_) {
            for (auto& bd : app.bounds) {
                if (bd)delete bd;
            }
        }
        remote_apps_.clear();

        // send and receive apps.
        int nerr = !intra_comm_ || intra_comm_->rank() == intra_root_
            ? send_recv_apps_()
            : 0;
        nerr = sync_intra_error_(nerr);
        if (nerr)return false;

        // 
        sync_field_info_();

        // allocate fields for coupled boundaries
        for (auto bd : data_.bounds) {
            bd->remove_all_field();
            for (auto& fd : data_.field_defs)
                bd->register_field(fd);
        }

        // allocate local fields for parallel solver of this application
        if (intra_comm_ && intra_comm_->size() > 1) {
            for (auto& bd : local_bounds_) {
                bd.remove_all_field();
                for (auto& fd : data_.field_defs)
                    bd.register_field(fd);
            }
        }
        
        // create interpolator
        create_interps_();

        is_started_ = true;
        return true;
    }

    void Application::stop_coupling()
    {
        if (!is_started_)return;

        if (!intra_comm_ || intra_comm_->rank() == intra_root_) {
            if (inter_comm_) {
                info("\nstop coupling!\n");
                inter_comm_->disconnect();
            }
        }

        is_started_ = false;
    }

    void Application::send_outgoing_fields_(double time)
    {
        // Operations on root solver of this application
        if (!intra_comm_ || intra_comm_->rank() == intra_root_) {
            // send outgoing fields
            for (int remote = 0; remote < inter_comm_->size(); ++remote) {
                if (remote == data_.rank)continue;

                // send time
                inter_comm_->send(&time, 1, remote, 0);

                // send fields
                auto& remote_app = remote_apps_.at(remote);
                for (auto& remote_fd : remote_app.field_defs) {
                    if (remote_fd.remote_app_rank != data_.rank)continue;

                    auto& this_fd = data_.field_defs.at(remote_fd.remote_field);
                    info("  sending outgoing field \"%s\" to %s\n", this_fd.name.c_str(), remote_app.app_name.c_str());

                    for (auto bd : data_.bounds) {
                        auto& f = bd->fields_.at(remote_fd.remote_field);
                        ASSERT(f.info == &this_fd);
                        inter_comm_->send(f.data.data(), static_cast<int>(f.data.numel()), remote, f.info->id);
                    }
                }
            }
        }
    }

    void Application::recv_incoming_fields_([[maybe_unused]]double time)
    {
        // Operations on root solver of this application
        if (!intra_comm_ || intra_comm_->rank() == intra_root_) {

            // receiving incoming fields
            for (int i = 0; i < inter_comm_->size(); ++i) {
                if (i == data_.rank)continue;

                auto& remote_app = remote_apps_.at(i);

                // recv time
                inter_comm_->recv(&remote_app.time, 1, i, 0);

                // recv fields
                for (auto& this_fd : data_.field_defs) {
                    if (this_fd.remote_app_rank != i)continue;

                    info("  receiving incoming field \"%s\" from %s\n", this_fd.name.c_str(), remote_app.app_name.c_str());

                    for (auto bd : remote_app.bounds) {
                        auto& f = bd->fields_.at(this_fd.remote_field);
                        inter_comm_->recv(f.data.data(), static_cast<int>(f.data.numel()), i, f.info->id);
                    }
                }
            }

            // TODO: time interpolate

            // interpolate incoming fields
            this_fields_.resize(data_.bounds.size(), nullptr);
            for (auto& this_fd : data_.field_defs) {
                if (this_fd.is_orphan || (this_fd.iotype != IncomingDofs && this_fd.iotype != IncomingLoads))continue;

                // find interpolator
                Interpolator* ip = field_interps_.at(this_fd.id);
                if (!ip) {
                    error("failed get interpolator!!!");
                    return;
                }

                // this fields
                for (size_t i = 0; i < data_.bounds.size(); ++i) {
                    auto bd = data_.bounds.at(i);
                    this_fields_.at(i) = &bd->fields_.at(this_fd.id);
                }
                //remote fields
                auto& remote_app = remote_apps_.at(this_fd.remote_app_rank);
                remote_fields_.resize(remote_app.bounds.size(), nullptr);
                for (size_t i = 0; i < remote_app.bounds.size(); ++i) {
                    auto bd = remote_app.bounds.at(i);
                    remote_fields_.at(i) = &bd->fields_.at(this_fd.remote_field);
                }

                // interpolate DOFs from remote
                if (ip->source_app_id() == remote_app.rank && this_fd.iotype == IncomingDofs) {
                    ip->interp_dofs_s2t(
                        make_span(remote_fields_.data(), remote_fields_.size()),
                        make_span(this_fields_.data(), this_fields_.size())
                    );
                }
                // interpolate LOADs for this
                else if (this_fd.iotype == IncomingLoads) {
                    ip->interp_load_t2s(
                        make_span(remote_fields_.data(), remote_fields_.size()),
                        make_span(this_fields_.data(), this_fields_.size())
                    );
                }
                else {
                    error("invalid interpolator!!!");
                    return;
                }
            }
        }
    }

    void Application::exchange_solution(double time, void* user_data)
    {
        if (data_.time != time) { data_.iter = 0; }
        data_.time = time; ++data_.iter;

        read_outgoing_(user_data);
        send_outgoing_fields_(time);
        recv_incoming_fields_(time);
        write_incoming_(user_data);
    }
    
    int Application::send_recv_apps_()
    {
        if (intra_comm_ && intra_comm_->rank() != intra_root_)return 0;

        data_.rank = inter_comm_->rank();

        int nerr = 0;

        //--- Phase 1: sending application to other without boundary data

        for (int i = 0; i < inter_comm_->size(); ++i) {
            if (i == inter_comm_->rank())continue;

            int tag = 0;
            int buf[2] = {
                static_cast<int>(data_.bounds.size()),
                static_cast<int>(data_.field_defs.size())
            };

            // application name
            inter_comm_->send(data_.app_name, i, tag);

            // boundary and field number
            inter_comm_->send(buf, 2, i, ++tag);
        }

        // receive application
        int n_remote_bd = 0;
        remote_apps_.resize(inter_comm_->size());
        for (int i = 0; i < inter_comm_->size(); ++i) {
            if (i == inter_comm_->rank())continue;

            auto& app = remote_apps_.at(i);
            app.rank  = i;

            int tag = 0;
            int buf[2] = { 0 };

            // application name
            inter_comm_->recv(app.app_name, i, tag);

            // boundary and field number
            inter_comm_->recv(buf, 5, i, ++tag);

            // allocate
            app.bounds.resize(buf[0], nullptr);
            app.field_defs.resize(buf[1]);

            n_remote_bd += buf[0];
        }

        // Phase 2: allocate remote boundaries.
        remote_bounds_.resize(n_remote_bd);

        // Phase 3: send and receive boundary and field info.

        for (int i = 0; i < inter_comm_->size(); ++i) {
            if (i == inter_comm_->rank())continue;

            int tag = 0;

            // boundary
            for (Boundary* bd : data_.bounds)inter_comm_->send(*bd, i, ++tag);

            // field info
            for (auto& info : data_.field_defs)inter_comm_->send(info, i, ++tag);
        }
        for (int i = 0, ib = 0; i < inter_comm_->size(); ++i) {
            if (i == inter_comm_->rank())continue;

            auto& app = remote_apps_.at(i);

            int tag = 0;

            // boundary
            for (auto& bd : app.bounds) {
                bd = &remote_bounds_.at(ib);
                inter_comm_->recv(*bd, i, ++tag);
                ++ib;
            }

            // fields
            for (auto& info : app.field_defs)
                inter_comm_->recv(info, i, ++tag);
        }

        //----------------------------------------------
        // Check fields
        //----------------------------------------------

        // reset all fields as orphan and out of date
        for (auto& app : remote_apps_) {
            for (auto & fd : app.field_defs) {
                fd.is_orphan       = true;
                fd.is_out_of_date  = true;
                fd.remote_app_rank = -1;
                fd.id              = static_cast<int>(&fd - app.field_defs.data());
            }
        }
        // find source and target fields for this application
        for (auto & f : data_.field_defs) {
            f.is_orphan       = true;
            f.is_out_of_date  = true;
            f.remote_app_rank = -1;
            f.id              = static_cast<int>(&f - data_.field_defs.data());
            for (auto& app : remote_apps_) {
                for (auto & rf : app.field_defs) {
                    if (f.name != rf.name)continue;
                    if ((f.iotype == IncomingDofs  && rf.iotype == OutgoingDofs) ||
                        (f.iotype == IncomingLoads && rf.iotype == OutgoingLoads) ||
                        (f.iotype == OutgoingDofs  && rf.iotype == IncomingDofs) ||
                        (f.iotype == OutgoingLoads && rf.iotype == IncomingLoads)) {
                        f.is_orphan        = false;
                        f.remote_app_rank  = app.rank;
                        f.remote_field     = rf.id;
                        rf.is_orphan       = false;
                        rf.remote_app_rank = data_.rank;
                        rf.remote_field    = f.id;
                        break;
                    }
                }
            }
            if (f.is_orphan && (f.iotype == IncomingDofs || f.iotype == IncomingLoads)) {
                info("\n***ERROR*** source for incoming filed \"%s\" is not found!\n", f.name.c_str());
                ++nerr;
            }
            if (f.is_orphan && (f.iotype == OutgoingDofs || f.iotype == OutgoingLoads)) {
                info("\n***WARNING*** target for outgoing filed \"%s\" is not found!\n", f.name.c_str());
            }
        }

        // sync error for all applications
        for (int i = 0; i < inter_comm_->size(); ++i) {
            if (i == inter_comm_->rank())continue;
            inter_comm_->send(&nerr, 1, i, 0);
        }
        for (int i = 0; i < inter_comm_->size(); ++i) {
            if (i == inter_comm_->rank())continue;
            int n = 0;
            inter_comm_->recv(&n, 1, i, 0);
            nerr += n;
        }

        return nerr;
    }

    int Application::sync_intra_error_(int nerr)
    {
        // send and 
        if (intra_comm_ && intra_comm_->size() > 0) {
            if (intra_comm_->rank() == intra_root_) {
                for (int i = 0; i < intra_comm_->size(); ++i) {
                    if (i != intra_root_)
                        intra_comm_->send(&nerr, 1, i, 0);
                }
            }
            else {
                intra_comm_->recv(&nerr, 1, intra_root_, 0);
            }
        }
        return nerr;
    }

    void Application::sync_field_info_()
    {
        if (!intra_comm_ || intra_comm_->size() <= 1)return;

        struct Info{
            int
                remote_app_rank,
                remote_field,
                orphan;
        };
        int nf = static_cast<int>(data_.field_defs.size());
        std::vector<Info> data(nf);
        if (intra_comm_->rank() == intra_root_) {
            for (size_t i = 0; i < data.size(); ++i) {
                auto& d = data.at(i);
                auto& fd = data_.field_defs.at(i);
                d.remote_app_rank = fd.remote_app_rank;
                d.remote_field = fd.remote_field;
                d.orphan = fd.is_orphan;
            }
            for (int i = 0; i < intra_comm_->size(); ++i) {
                if (i == intra_comm_->rank())continue;
                intra_comm_->send((const int*)data.data(), 3 * nf, i, 0);
            }
        }
        else {
            intra_comm_->recv((int*)data.data(), 3 * nf, intra_root_, 0);

            for (size_t i = 0; i < data.size(); ++i) {
                auto& d = data.at(i);
                auto& fd = data_.field_defs.at(i);
                fd.remote_app_rank = d.remote_app_rank;
                fd.remote_field = d.remote_field;
                fd.is_orphan = d.orphan;
            }
        }
    }

    void Application::create_interps_()
    {
        if (intra_comm_ && intra_comm_->rank() != intra_root_)return;

        interps_.clear();
        interps_.reserve(data_.field_defs.size());
        field_interps_.resize(data_.field_defs.size());

        for (auto& fd : data_.field_defs) {
            // skip orphan field
            if (fd.is_orphan)continue;

            //? this application is used as source if contains any outgoing DOFs.
            if (fd.iotype == OutgoingDofs) {
                int src = data_.rank;
                int des = fd.remote_app_rank;
                for (auto& ip : interps_) {
                    if (ip.source_app_id() == src && ip.target_app_id() == des) {
                        field_interps_.at(fd.id) = &ip;
                        break;
                    }
                }
                if (field_interps_.at(fd.id))continue;

#if __cplusplus >= 201703L
                auto& ip = interps_.emplace_back(Interpolator{});
#else
                interps_.push_back(Interpolator{});
                auto& ip = interps_.back();
#endif
                ip.set_app_id(src, des);

                for (auto bd : data_.bounds)ip.add_source_boundary(*bd);
                for (auto bd : remote_apps_.at(fd.remote_app_rank).bounds)
                    ip.add_target_boundary(*bd);

                ip.compute_interp_coeff();
            }
            //? We do not create interpolator if application does contain any outgoing DOFs.
            // remote application is used as source
            else {
                int des = data_.rank;
                int src = fd.remote_app_rank;
                for (auto& ip : interps_) {
                    if (ip.source_app_id() == src && ip.target_app_id() == des) {
                        field_interps_.at(fd.id) = &ip;
                        break;
                    }
                }
                if (field_interps_.at(fd.id))continue;

#if __cplusplus >= 201703L
                auto& ip = interps_.emplace_back(Interpolator{});
#else
                interps_.push_back(Interpolator{});
                auto& ip = interps_.back();
#endif
                ip.set_app_id(src, des);

                for (auto bd : remote_apps_.at(fd.remote_app_rank).bounds)
                    ip.add_source_boundary(*bd);
                for (auto bd : data_.bounds)ip.add_target_boundary(*bd);

                ip.compute_interp_coeff();
            }
        }
    }

    void Application::read_outgoing_(void* user_data)
    {
        // Parallel running: read and assemble outgoing fields
        if (intra_comm_ && intra_comm_->size() > 1) {
            // reading data to local-field
            for (size_t i = 0; i < data_.field_defs.size(); ++i) {
                auto& fd = data_.field_defs.at(i);
                if (fd.is_orphan || (fd.iotype != OutgoingDofs && fd.iotype != OutgoingLoads))continue;

                if (intra_comm_->rank() == intra_root_)
                    info("  reading outgoing field: %s\n", fd.name.c_str());

                for (size_t ib = 0; ib < local_bounds_.size(); ++ib) {
                    auto& local_bd = local_bounds_.at(ib);
                    auto& local_f  = local_bd.fields_.at(i);
                    ASSERT(local_f.info == &fd);
                    getter_(this, &local_bd, fd.name.c_str(), fd.ncomp, fd.location, local_f.data.data(), user_data);

                    // assemble the full-field
                    auto& global_bd = bounds_.at(ib).full_boundary();
                    auto& global_f  = global_bd.fields_.at(i);
                    ASSERT(global_f.info == &fd);

                    if (fd.location == NodeCentered)
                        bounds_.at(ib).gather_node_fields(fd.ncomp, local_f.data.data(), global_f.data.data());
                    else
                        bounds_.at(ib).gather_face_fields(fd.ncomp, local_f.data.data(), global_f.data.data());
                }                
            }
        }
        // Serial running: read outgoing fields
        else {
            // reading data to full-field
            for (size_t i = 0; i < data_.field_defs.size(); ++i) {
                auto& fd = data_.field_defs.at(i);
                if (fd.is_orphan || (fd.iotype != OutgoingDofs && fd.iotype != OutgoingLoads))continue;

                info("  reading outgoing field: %s\n", fd.name.c_str());

                for (auto bd : data_.bounds) {
                    auto& f = bd->fields_.at(i);
                    ASSERT(f.info == &fd);
                    getter_(this, bd, fd.name.c_str(), fd.ncomp, fd.location, f.data.data(), user_data);
                }
            }
        }
    }
    void Application::write_incoming_(void* user_data)
    {
        // Parallel running: scatter and writing incoming fields
        if (intra_comm_ && intra_comm_->size() > 1) {
            for (size_t i = 0; i < data_.field_defs.size(); ++i) {
                auto& fd = data_.field_defs.at(i);
                if (fd.is_orphan || (fd.iotype != IncomingDofs && fd.iotype != IncomingLoads))continue;

                if (intra_comm_->rank() == intra_root_)
                    info("  writing incoming field: %s\n", fd.name.c_str());

                for (size_t ib = 0; ib < local_bounds_.size(); ++ib) {
                    auto& local_bd = local_bounds_.at(ib);
                    auto& local_f = local_bd.fields_.at(i);
                    ASSERT(local_f.info == &fd);
                    auto& global_bd = bounds_.at(ib).full_boundary();
                    auto& global_f = global_bd.fields_.at(i);
                    ASSERT(global_f.info == &fd);

                    // scatter fields
                    if (fd.location == NodeCentered)
                        bounds_.at(ib).scatter_node_fields(fd.ncomp, global_f.data.data(), local_f.data.data());
                    else
                        bounds_.at(ib).scatter_face_fields(fd.ncomp, global_f.data.data(), local_f.data.data());

                    // writing
                    setter_(this, &local_bd, fd.name.c_str(), fd.ncomp, fd.location, local_f.data.data(), user_data);
                }
            }
        }
        // Serial running: writing incoming fields
        else {
            for (size_t i = 0; i < data_.field_defs.size(); ++i) {
                auto& fd = data_.field_defs.at(i);
                if (fd.is_orphan || (fd.iotype != IncomingDofs && fd.iotype != IncomingLoads))continue;

                if (!intra_comm_ || intra_comm_->rank() == intra_root_)
                    info("  writing incoming field: %s\n", fd.name.c_str());

                for (auto bd : data_.bounds) {
                    auto& f = bd->fields_.at(i);
                    ASSERT(f.info == &fd);
                    setter_(this, bd, fd.name.c_str(), fd.ncomp, fd.location, f.data.data(), user_data);
                }
            }
        }
    }

    void Application::save_tecplot(const char* file, bool without_fields/* = false*/)
    {
        if (intra_comm_ && intra_comm_->rank() != intra_root_)return;

        std::ofstream ofs(file);
        if (!ofs.is_open()) {
            warn("failed opening Tecplot file: %s", file);
            return;
        }

        // TITLE
        ofs << "TITLE = \"" << data_.app_name << "\"\n";

        // FILETYPE
        ofs << "FILETYPE = \"" << (without_fields ? "GRID" : "FULL") << "\"\n";

        // VARIABLES
        ofs << "\"X\" \"Y\" \"Z\"";
        std::string var_loc;
        int nvar = 3;
        if (!without_fields) {
            for (auto& fd : data_.field_defs) {
                if (fd.is_orphan)continue;

                if (fd.ncomp == 1) {
                    ofs << " \"" << fd.name << "\"";
                    ++nvar;

                    // ,xx
                    if (var_loc.empty())
                        var_loc = "([";
                    else
                        var_loc.push_back(',');
                    var_loc.append(std::to_string(nvar));
                }
                else {
                    for (int j = 0; j < fd.ncomp; ++j) {
                        ofs << " \"" << fd.name << '[' << j + 1 << "]\"";
                        ++nvar;

                        // ,xx
                        if (var_loc.empty())
                            var_loc = "([";
                        else
                            var_loc.push_back(',');
                        var_loc.append(std::to_string(nvar));
                    }
                }
            }
            if (!var_loc.empty())var_loc.append("]=CELLCENTERED)");
        }
        ofs << '\n';

        auto write_field = [](std::ofstream& ofs, const double* data, int count, int stride = 1) {
            int i = 1;
            for (; i <= count; ++i) {
                ofs << ' ' << std::scientific << *data;
                if (i % 8 == 0)ofs << '\n';
                data += stride;
            }
            if (i % 8 != 0)ofs << '\n';
        };

        // ZONES
        for (auto bd : data_.bounds) {
            // zone title
            ofs << "ZONE\n"
                << "T = \"" << bd->name() << "\"\n";

            // 
            bool all_tri  = bd->face_type_num(TRI3) + bd->face_type_num(TRI6) == bd->nface();
            //bool all_quad = bd->face_count(FT_QUAD4) + bd->face_count(FT_QUAD8) == bd->nface();

            // Polygon not exists: write finite element grid (ignore middle node)
            if (!bd->contains_polygon()) {
                switch (bd->topo()) {
                case ZT_POINTS:  ofs << "ZONETYPE = ORDERED\n"         << "I = "     << bd->nnode() << '\n'; break;
                case ZT_CURVE:   ofs << "ZONETYPE = FELINESEG\n"       << "NODES = " << bd->nnode() << " ELEMENTS = " << bd->nface() << '\n'; break;
                case ZT_SURFACE:
                    if (all_tri )
                        ofs << "ZONETYPE = FETRIANGLE\n" << "NODES = " << bd->nnode() << " ELEMENTS = " << bd->nface() << '\n';
                    else
                        ofs << "ZONETYPE = FEQUADRILATERAL\n" << "NODES = " << bd->nnode() << " ELEMENTS = " << bd->nface() << '\n';
                    break;
                //case ZT_VOLUME:  ofs << "ZONETYPE = FEBRICK\n"         << "NODES = " << bd->nnode() << " ELEMENTS = " << bd->nface() << '\n'; break;
                }
            }
            // Polygon exists: write face-based grid
            else {
                bd->create_edges(); // create edges

                ofs << "ZONETYPE = FELINESEG\n"
                    << "NODES = " << bd->nnode()
                    << " FACES = " << bd->edges_for_surface().size()
                    << " ELEMENTS = " << bd->nface() << '\n';
            }

            // DATATYPE
            ofs << "DT = DOUBLE DATAPACKING = BLOCK\n";

            // VARLOCATION
            if (!var_loc.empty())ofs << "VARLOCATION = " << var_loc << '\n';

            // DATA
            write_field(ofs, bd->node_coords().data()->data() + 0, bd->nnode(), 3);
            write_field(ofs, bd->node_coords().data()->data() + 1, bd->nnode(), 3);
            write_field(ofs, bd->node_coords().data()->data() + 2, bd->nnode(), 3);

            // fields
            for (auto& fd : data_.field_defs) {
                if (fd.is_orphan)continue;

                auto& field = bd->fields().at(fd.id);
                for (int j = 0; j < fd.ncomp; ++j) {
                    write_field(ofs, field.data.data() + j, (int_l)field.data.nrow(), fd.ncomp);
                }
            }

            // element data is ignored for points zone.
            if (bd->topo() == ZT_POINTS)continue;

            // face-nodes connections
            // FELINESEG zone
            if (bd->topo() == ZT_CURVE) {
                for (int_l i = 0; i < bd->nface(); ++i) {
                    auto&& nodes = bd->face_nodes()[i];
                    ofs << nodes[0] + 1 << ' ' << nodes[1] + 1 << '\n'; //? convert to one-based index, ignore middle node for FT_BAR3
                }
            }
            // triangle or quadrilateral
            else if (!bd->contains_polygon()) {
                if (all_tri) {
                    for (int_l i = 0; i < bd->nface(); ++i) {
                        auto&& nodes = bd->face_nodes()[i];
                        ofs << nodes[0] + 1 << ' ' << nodes[1] + 1 << ' ' << nodes[2] + 1 << '\n'; //? convert to one-based index, ignore middle node for FT_TRI6
                    }
                }
                else {
                    for (int_l i = 0; i < bd->nface(); ++i) {
                        auto&& nodes = bd->face_nodes()[i];
                        switch (bd->face_types().at(i)) {
                        case TRI3: ofs << nodes[0] + 1 << ' ' << nodes[1] + 1 << ' ' << nodes[2] + 1 << '\n'; break;
                        case TRI6: ofs << nodes[0] + 1 << ' ' << nodes[1] + 1 << ' ' << nodes[2] + 1 << '\n'; break;
                        case QUAD4:ofs << nodes[0] + 1 << ' ' << nodes[1] + 1 << ' ' << nodes[2] + 1 << ' ' << nodes[3] + 1 << '\n'; break;
                        case QUAD8:ofs << nodes[0] + 1 << ' ' << nodes[1] + 1 << ' ' << nodes[2] + 1 << ' ' << nodes[3] + 1 << '\n'; break;
                        }
                    }
                }
            }
            // polygon surface
            else {
                // face nodes
                for (auto& e : bd->edges_for_surface()) {
                    ofs << e.n0 + 1 << ' ' << e.n1 + 1 << '\n';
                }
                // left element
                for (auto& e : bd->edges_for_surface()) {
                    auto f0 = e.f0 != invalid_id ? e.f0 + 1 : 0;
                    ofs << ' ' << f0;
                }
                ofs << '\n';
                // right element
                for (auto& e : bd->edges_for_surface()) {
                    auto f1 = e.f1 != invalid_id ? e.f1 + 1 : 0;
                    ofs << ' ' << f1;
                }
                ofs << '\n';
            }
        }

        ofs.close();
    }
}
