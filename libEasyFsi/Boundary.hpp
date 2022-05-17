#pragma once
#include <istream>
#include <ostream>
#include <vector>
#include <string>
#include <array>

#include "Index.hpp"
#include "TinyVector.hpp"
#include "TinyMatrix.hpp"
#include "IndexSet.hpp"
#include "MeshConnectivity.hpp"
#include "KDTree.hpp"
#include "DynamicMatrix.hpp"
#include "Field.hpp"

namespace EasyLib {

    //! @brief Topology type of face on coupled boundary.
    typedef enum FaceTopo
    {
        FT_BAR2 = 0,//! 2-nodes linear element
        FT_BAR3,    //! 3-nodes quadratic element
        FT_TRI3,    //! 3-nodes linear triangle element
        FT_TRI6,    //! 6-nodes quadratic triangle element
        FT_QUAD4,   //! 4-node linear quadrilateral element
        FT_QUAD8,   //! 8-node quadratic quadrilateral element
        FT_POLYGON  //! general polygon element (node number > 4)
    }FaceTopo;

    //! @brief Node number of specific face topology.
    extern const int npf[FT_POLYGON + 1];

    //! @brief Maximum number of face.
    inline constexpr int npf_max = 8;

    //! @brief The order of face, 1=Linear£¬2=Quadratic
    extern const int face_order[FT_POLYGON + 1];

    //! @brief Topology of coupled zone.
    typedef enum ZoneTopo
    {
        ZT_POINTS = 0, //! Discrete points
        ZT_CURVE,      //! Curve Zone
        ZT_SURFACE,    //! Surface Zone
        ZT_VOLUME      //! Volume Zone
    }ZoneTopo;

    //! @brief Shape of coupled zone.
    typedef enum ZoneShape
    {
        ZS_POINT = 0, //! All nodes are coincide.
        ZS_COLINEAR,  //! All nodes are colinear.
        ZS_COPLANER,  //! All nodes are coplaner.
        ZS_GENERAL    //! General 3D distribution.
    }ZoneShape;

    //! @brief Coupled boundary class.
    class Boundary
    {
    public:
        using dvec = std::vector<double>;
        using ivec = std::vector<int>;
        using vec3 = TinyVector<double, 3>;
        using vvec = std::vector<vec3>;

        struct Edge
        {
            int_l n0{ 0 }, n1{ 0 }, f0{ 0 }, f1{ 0 };

            inline friend bool operator < (const Edge& lhs, const Edge& rhs)
            {
                return lhs.n0 < rhs.n0 || (lhs.n0 == rhs.n0 && lhs.n1 < rhs.n1);
            }
        };

        Boundary() = default;
        Boundary(const Boundary&) = default;
        Boundary& operator = (const Boundary&) = default;
        Boundary(Boundary&& bd)noexcept;
        Boundary& operator = (Boundary&& bd)noexcept;

        //! @brief Clear all boundary data.
        void clear();

        //! @brief Preallocate memory for fast insertion.
        //! @param max_node   Maximum node number
        //! @param max_face   Maximum face number
        //! @param max_fnodes Maximum face-nodes number
        void reserve(int_l max_node, int_l max_face, int_l max_fnodes);

        void set_name(const char* sname);

        //! @brief Add a node to boundary, and return it's local index.
        //! @param x Coordinate-x of the node
        //! @param y Coordinate-y of the node
        //! @param z Coordinate-z of the node
        //! @param global_id Global unique index of the node, used to distinguish node.
        //! @note Return old id if node with same global index already exists.
        int add_node(double x, double y, double z, int_g global_id = -1);

        //! @brief Add a node to boundary, and return it's local index.
        //! @param coord      Coordinates of the node
        //! @param global_id  Global unique index of the node, used to distinguish node.
        //! @note Return old id if node with same global index already exists.
        int add_node(const vec3& coord, int_g global_id = -1);

        //! @brief Add a face to boundary, and return it's local index.
        //! @param type    Face shape.
        //! @param count   Node number of face, should be \npf[\type].
        //! @param fnodes  Local indices of face nodes.
        //! @return Return face index.
        int add_face(FaceTopo type, int count, const int_l* fnodes);
        int add_face(FaceTopo type, int count, const int_l* fnodes, double cx, double cy, double cz);
        int add_face(FaceTopo type, int count, const int_l* fnodes, const vec3& fcent);

        //! @brief Modify node coordinates.
        //! @param id  Index of the node
        //! @param x   New x-coordinate of the node.
        //! @param y   New y-coordinate of the node.
        //! @param z   New z-coordinate of the node.
        void set_node_coords(int_l id, double x, double y, double z);

        //! @brief Modify face centroid. Used to keep consistent with solver.
        //! @param face Index of the face
        //! @param cx   New x-coordinate of the centroid.
        //! @param cy   New y-coordinate of the centroid.
        //! @param cz   New z-coordinate of the centroid.
        void set_face_cent(int_l face, double cx, double cy, double cz);

        //! @brief Modify face area. Used to keep consistent with solver.
        //! @param face Index of the face
        //! @param sx   New x-component of the area vector.
        //! @param sy   New y-component of the area vector.
        //! @param sz   New z-component of the area vector.
        void set_face_area(int_l face, double sx, double sy, double sz);

        //! @brief Create edges for surface boundary.
        void create_edges(/*std::vector<std::pair<int_l, int_l> >& edge_nodes, std::vector<std::pair<int_l, int_l> >& edge_faces*/);

        //! @brief Compute area, normal, centroid of faces, and topo type.
        //! @note This member function should be called after add all nodes and faces, or after modifying all node coordinates.
        void compute_metics(double biased_angle_deg = 5);

        //! @brief Get kdtree of nodes. Used to search nearest nodes.
        //! @note KDTree is created by \compute_metrics, so it is available only after calling \compute_metrics.
        inline auto& kdtree()const { ASSERT(!mesh_changed_); return kdtree_; }

        //! @brief Compute global IPS/TPS/SPLINE matrix.
        void compute_global_xps_matrix();

        //! @brief Compute global IPS/TPS/SPLINE interpolation coefficients for one point.
        //! @param [in]  p      Coordinates of the target point.
        //! @param [out] coeff  Coefficients used to do interpolating.
        void compute_global_xps_interp_coeff(const vec3& p, std::span<double> coeff);

        //! @brief Compute local IPS/TPS/SPLINE interpolation coefficients for one point.
        //! @param [in]  p           Coordinates of the target point.
        //! @param [in]  max_donor   Maximum number of donor points
        //! @param [out] ids         Index list of donor nodes.
        //! @param [out] coeff       Coefficients used to do interpolating.
        //! @param [out] n_donor     Actual number of donor points.
        //! @param [in]  min_dist_sq The minimum squared-distance to determine whether two points coincide.
        void compute_local_xps_interp_coeff(const vec3& p, int max_donor, std::span<int_l> ids, std::span<double> coeff, int& n_donor, double min_dist_sq = 1E-20);

        //! @brief Compute interpolation coefficients for one point using element projection method.
        //! @param [in]  p           Coordinates of the target point.
        //! @param [out] ids         Index list of donor nodes.
        //! @param [out] coeff       Coefficients used to do interpolating.
        //! @param [out] n_donor     Number of donor points.
        //! @param [in]  min_dist_sq The minimum squared-distance to determine whether two points coincide.
        void compute_project_interp_coeff(const vec3& p, int_l(&ids)[npf_max], double(&coeff)[npf_max], int& n_donor, double min_dist_sq = 1E-20);

        //void save(const char* file, const double* node_disp = nullptr)const;

        //void read(const char* file);

        //void read_init_coord(const char* file);

        void read_gmsh(const char* file);

        inline auto& name()const noexcept { return name_; }

        inline auto topo ()const { ASSERT(!mesh_changed_); return topo_; }
        inline auto shape()const { ASSERT(!mesh_changed_); return shape_; }

        inline bool contains_polygon()const noexcept { return face_count_[FT_POLYGON] > 0; }
        inline bool contains_high_order_face()const noexcept { return face_count_[FT_BAR3] || face_count_[FT_TRI6] || face_count_[FT_QUAD8]; }
        //inline bool is_all_high_order()const { return is_high_order_; }

        inline auto& edges_for_surface()const noexcept { return edges_; }

        inline const vec3& coords_min()const noexcept { ASSERT(!mesh_changed_); return coord_min_; }
        inline const vec3& coords_max()const noexcept { ASSERT(!mesh_changed_); return coord_max_; }

        inline int_l node_num()const { return nodes_.size(); }
        inline int_l face_num()const { return face_nodes_.nrow(); }

        inline const IndexSet&         nodes     ()const noexcept { return nodes_; }
        inline const MeshConnectivity& face_nodes()const noexcept { return face_nodes_; }

        inline const auto& node_coords   ()const noexcept { return node_coords_; }
        inline const auto& face_centroids()const noexcept { ASSERT(!mesh_changed_); return face_centroids_; }
        inline const auto& face_areas    ()const noexcept { ASSERT(!mesh_changed_); return face_area_; }
        inline const auto& face_normals  ()const noexcept { ASSERT(!mesh_changed_); return face_normal_; }
        inline const auto& face_types    ()const noexcept { return face_types_; }

        inline auto& face_count(           )const noexcept { return face_count_; }
        inline auto  face_count(FaceTopo ft)const noexcept { return face_count_[ft]; }

        inline auto& get_fields()const noexcept { return fields_; }

        Field&       get_field(const char* field_name);
        const Field& get_field(const char* field_name)const;

        friend std::istream& operator >>(std::istream& is,       Boundary& bd);
        friend std::ostream& operator <<(std::ostream& os, const Boundary& bd);

        friend class Communicator;
        friend class Interpolator;
        friend class Application;

    private:
        void register_field_(const FieldInfo& fd);
        
    private:
        std::string      name_;
        IndexSet         nodes_;
        MeshConnectivity face_nodes_, node_faces_;
        vvec             node_coords_;
        vvec             face_centroids_;
        ivec             face_types_;
        dvec             face_area_;
        vvec             face_normal_;
        std::array<int_l, FT_POLYGON + 1> face_count_{ 0 };

        //---------------------------------------------------
        // Following Data will be ignored by communicator.
        //---------------------------------------------------

        bool             mesh_changed_{ false };//! Is mesh changed?

        ZoneTopo         topo_ { ZT_POINTS };   //! update in add_face()
        ZoneShape        shape_{ ZS_GENERAL };  //! update in compute_metics()

        std::vector<Edge> edges_;  //! Edges data for polygon faces.

        vec3              coord_min_{ 0,0,0 };  //! update in compute_metics()
        vec3              coord_max_{ 0,0,0 };  //! update in compute_metics()

        Fields                   fields_; //! for application

        KDTree<double, 3, int_l> kdtree_; //! 

        //--- data members used to compute global IPS/TPS/SPLINE matrix

        TinyMatrix<double, 4> xps_tm_;     //! used to transform point to local CS. update in compute_metics()
        dvec                  xps_coords_; //! 
        DynamicMatrix         xps_ts_inv_; //!
        bool                  xps_computed_{ false };

        //--- data members used to compute local IPS/TPS/SPLINE matrix

        vvec                  local_xps_points_;

        //--- buffers used by global/local IPS/TPS/SPLINE methods

        std::vector<int>      ibuffer_; //! used for inverse xps_ts_inv_.
        dvec                  dbuffer_; //! used for computing XPS interpolation coefficient.
    };
}
