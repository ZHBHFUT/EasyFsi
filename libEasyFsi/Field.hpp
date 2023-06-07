#pragma once
#include <string>

#include "DynamicMatrix.hpp"
#include "Logger.hpp"

namespace EasyLib {

    enum FieldLocation
    {
        NodeCentered = 0,
        FaceCentered,
        CellCentered
    };
    enum FieldIO
    {
        IncomingDofs = 0,
        IncomingLoads,
        OutgoingDofs,
        OutgoingLoads
    };

    struct FieldInfo
    {
        std::string   name;
        std::string   units;
        int           ncomp{ 0 };
        FieldLocation location{ NodeCentered };
        FieldIO       iotype{ IncomingDofs };

        int           id{ 0 };
        int           remote_app_rank{ -1 };
        int           remote_field{ -1 };
        double        time{ 0 };               // Timestamp for current values.
        bool          is_orphan{ false };      // This field is not used by any remote application.
        bool          is_out_of_date{ false }; // Field is out of date.

        inline bool is_incoming()const noexcept { return iotype == IncomingDofs || iotype == IncomingLoads; }
        inline bool is_outgoing()const noexcept { return iotype == OutgoingDofs || iotype == OutgoingLoads; }
        inline bool is_node_centered()const noexcept { return location == NodeCentered; }
    };

    struct Field
    {
        const FieldInfo* info{ nullptr };
        DynamicMatrix    data;

        inline void set(int_l i, int_l j, double value)
        {
            if (info && !info->is_outgoing())
                error("unable write field because it is not outgoing field!");
            data.at(i, j) = value;
        }
        inline double get(int_l i, int_l j)const
        {
            return data.at(i, j);
        }
    };

    using Fields = std::vector<Field>;
}
