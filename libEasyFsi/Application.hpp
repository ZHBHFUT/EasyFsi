#pragma once
#include <string>

#include "Boundary.hpp"
#include "Communicator.hpp"
#include "DistributedBoundary.hpp"
#include "Interpolator.hpp"

namespace EasyLib {

    typedef void(__stdcall *get_boundary_field_function)(const Boundary* bd, const char* name, int ncomp, FieldLocation loc, double*       data, void* user_data);
    typedef void(__stdcall *set_boundary_field_function)(const Boundary* bd, const char* name, int ncomp, FieldLocation loc, const double* data, void* user_data);

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

        //void recv_incoming();
        //void send_outgoing();

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

    private:
        struct ApplicationData
        {
            std::string            app_name;
            DynamicArray<Boundary*, 1> bounds;
            int                    rank{ -1 };
            int                    iter{ 0 };
            double                 time{ 0 };
            DynamicArray<FieldInfo, 1> field_defs;
        };

        Communicator*                    inter_comm_{ nullptr }; //! communicator between coupled applications
        Communicator*                    intra_comm_{ nullptr }; //! communicator between each partition of current application.
        int                              intra_root_{ 0 };
        ApplicationData                  data_;

        DynamicArray<ApplicationData,1>      remote_apps_; //! available only on root process.
        DynamicArray<Boundary, 1>            local_bounds_;
        DynamicArray<DistributedBoundary, 1> bounds_;
        DynamicArray<Boundary, 1>            remote_bounds_;

        bool                             is_started_{ false };
        DynamicArray<Interpolator, 1>    interps_;
        DynamicArray<Interpolator*, 1>   field_interps_;

        get_boundary_field_function      getter_{ nullptr };
        set_boundary_field_function      setter_{ nullptr };

        DynamicArray<Field*, 1> this_fields_;
        DynamicArray<Field*, 1> remote_fields_;
    };
}
