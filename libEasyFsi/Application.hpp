#pragma once
#include <string>

#include "Boundary.hpp"
#include "Communicator.hpp"
#include "DistributedBoundary.hpp"
#include "Interpolator.hpp"

namespace EasyLib {

    struct ApplicationData
    {
        std::string            app_name;
        std::vector<Boundary*> bounds;
        int                    rank{ -1 };
        int                    iter{ 0 };
        double                 time{ 0 };
        //std::vector<Fields>    bd_fields;
        std::vector<FieldInfo> field_defs;
    };

    typedef void(*get_boundary_field_function)(const Boundary* bd, const char* name, int ncomp, FieldLocation loc, double*       data);
    typedef void(*set_boundary_field_function)(const Boundary* bd, const char* name, int ncomp, FieldLocation loc, const double* data);

    class Application
    {
    public:
        using dv = DynamicVector;

        explicit Application(const char* name);
        Application(const char* name, Communicator& intra_comm, int root);

        virtual ~Application() = default;

        void clear();

        void create(const char* name);
        void create(const char* name, Communicator& intra_comm, int root);

        Boundary& add_coupled_boundary();

        void set_field_function(get_boundary_field_function getter, set_boundary_field_function setter);

        void register_field(const char* field_name, int ncomp, FieldLocation location, FieldIO iotype, const char* units);

        bool start_coupling(Communicator& comm);

        void exchange_solution();

        void stop_coupling();

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

        void read_outgoing_();
        void write_incoming_();

    private:
        Communicator*                    inter_comm_{ nullptr }; //! communicator between coupled applications
        Communicator*                    intra_comm_{ nullptr }; //! communicator between each partition of current application.
        int                              intra_root_{ 0 };
        ApplicationData                  data_;
        std::vector<ApplicationData>     remote_apps_; //! only available on root process
        std::vector<Boundary>            local_bounds_;
        //std::vector<Fields>              local_bd_fields_;
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
