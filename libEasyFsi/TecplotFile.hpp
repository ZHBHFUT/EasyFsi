#pragma once
#include <string>
#include <vector>

#include "DynamicVector.hpp"
#include "MeshConnectivity.hpp"

namespace EasyLib {

    class TecplotFile;

    enum TecZoneType
    {
        ORDERED,
        FELINESEG,
        FETRIANGLE,
        FEQUADRILATERAL,
        FETETRAHEDRON,
        FEBRICK,
        FEPOLYGON,
        FEPOLYHEDRAL
    };
    enum TecVarLocation
    {
        NODAL,
        CELLCENTERED
    };

    struct TecplotZone
    {
        std::string title;
        TecZoneType type{ ORDERED };

        int
            i{ 0 }, j{ 0 }, k{ 0 };

        int
            nodes         { 0 },
            elements      { 0 },
            faces         { 0 },
            num_face_nodes{ 0 },
            num_bd_faces  { 0 };

        bool is_block{ false };

        double
            solution_time{ 0 };

        //std::vector<int> face_neighbor_mode_;

        std::vector<TecVarLocation> var_loc;
        std::vector<int>            var_share_list;
        std::vector<DynamicVector>  vars;

        MeshConnectivity                 elem_nodes; //! zero-based index
        MeshConnectivity                 face_nodes; //! zero-based index
        std::vector<std::pair<int, int>> face_elems; //! zero-based index
    };

    class TecplotFile
    {
    public:
        enum FileType
        {
            FULL,
            GRID,
            SOLUTION
        };

        void read_ascii(const char* file);


    private:
        std::string              title_;
        FileType                 type_{ FULL };
        std::vector<std::string> variables_;
        std::vector<TecplotZone> zones_;
    };

}
