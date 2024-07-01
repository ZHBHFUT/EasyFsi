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
//! @file       Application.hpp
//!             The definition of Application class
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <string>

#include "Boundary.hpp"
#include "Communicator.hpp"
#include "DistributedBoundary.hpp"
#include "Interpolator.hpp"

namespace EasyLib {
    class Application;

    typedef void(__stdcall *get_boundary_field_function)(const Application* app, const Boundary* bd, const char* name, int ncomp, FieldLocation loc, double*       data, void* user_data);
    typedef void(__stdcall *set_boundary_field_function)(const Application* app, const Boundary* bd, const char* name, int ncomp, FieldLocation loc, const double* data, void* user_data);

    class Application
    {
    public:
        using dv = DynamicVector;

        Application() = default;
        explicit Application(const char* name);
        Application(const char* name, Communicator& intra_comm, int root);

        virtual ~Application() = default;

        void clear();

        //void create(const char* name);
        //void create(const char* name, Communicator& intra_comm, int root);

        Boundary& add_coupled_boundary();

        int boundary_num()const { return static_cast<int>(bounds_.size()); }

        Boundary* boundary(int ib) { return ib >= 0 && ib < local_bounds_.size() ? &local_bounds_[ib] : nullptr; }

        void set_field_function(get_boundary_field_function getter, set_boundary_field_function setter);

        void register_field(const char* field_name, int ncomp, FieldLocation location, FieldIO iotype, const char* units);

        bool start_coupling(Communicator& comm);

        void exchange_solution(double time, void* user_data);

        void stop_coupling();

        void save_tecplot(const char* file, bool without_fields = false);

    private:
        //! @brief Send and receive application information to/from other applications
        //! @return 0=Succeed, other = error
        int  send_recv_apps_();

        //! @brief Create interpolator.
        void create_interps_();

        //! @brief Sync error between solver processes of this application.
        int  sync_intra_error_(int nerr);

        //! @brief Send field information to other solver processes of this application.
        void sync_field_info_();

        void read_outgoing_(void* user_data);
        void write_incoming_(void* user_data);

        //! @brief Read and send outgoing fields to other applications.
        //! @param time Current physical time.
        void send_outgoing_fields_(double time);

        //! @brief Receive and update incoming fields from other applications.
        //! @param time Current physical time.
        void recv_incoming_fields_(double time);

        void interp_incoming_();

    private:
        struct ApplicationData
        {
            std::string            app_name;  // name of application
            std::vector<Boundary*> bounds;    // boundaries
            int                    rank{ -1 };// 
            int                    iter{ 0 }; //
            double                 time{ 0 }; //
            std::vector<FieldInfo> field_defs;//
        };

        Communicator*                    inter_comm_{ nullptr }; //! communicator between coupled applications
        Communicator*                    intra_comm_{ nullptr }; //! communicator between each partition of current application.
        int                              intra_root_{ 0 };
        ApplicationData                  data_;

        std::vector<ApplicationData>     remote_apps_; //! available only on root process.
        std::vector<Boundary>            local_bounds_;
        std::vector<DistributedBoundary> bounds_;
        std::vector<Boundary>            remote_bounds_;

        bool                             is_started_{ false };
        std::vector<Interpolator>        interps_;
        std::vector<Interpolator*>       field_interps_;

        get_boundary_field_function      getter_{ nullptr };
        set_boundary_field_function      setter_{ nullptr };

        std::vector<Field*> this_fields_;
        std::vector<Field*> remote_fields_;
    };
}
