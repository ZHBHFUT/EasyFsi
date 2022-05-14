#pragma once
#include <vector>
//#include <set>

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
        //! @param [in]  p          Coordinates of the target point.
        //! @param [in]  max_donor  Maximum number of donor points
        //! @param [out] ids        Index list of donor nodes.
        //! @param [out] coeff      Coefficients used to do interpolating.
        //! @param [out] n_donor    Actual number of donor points.
        void compute_local_xps_interp_coeff(const vec3& p, int max_donor, std::span<int_l> ids, std::span<double> coeff, int& n_donor);

        //! @brief Compute interpolation coefficients for one point using element projection method.
        //! @param [in]  p        Coordinates of the target point.
        //! @param [out] ids      Index list of donor nodes.
        //! @param [out] coeff    Coefficients used to do interpolating.
        //! @param [out] n_donor  Number of donor points.
        void compute_project_interp_coeff(const vec3& p, int_l(&ids)[npf_max], double(&coeff)[npf_max], int& n_donor);

        //void save(const char* file, const double* node_disp = nullptr)const;

        //void read(const char* file);

        //void read_init_coord(const char* file);

        //void read_gmsh(const char* file);

        inline const vec3& coords_min()const noexcept { ASSERT(!mesh_changed_); return coord_min_; }
        inline const vec3& coords_max()const noexcept { ASSERT(!mesh_changed_); return coord_max_; }

        inline int_l node_num()const { return nodes_.size(); }
        inline int_l face_num()const { return face_nodes_.nrow(); }

        inline const IndexSet&         nodes     ()const { return nodes_; }
        inline const MeshConnectivity& face_nodes()const { return face_nodes_; }

        inline const auto& node_coords   ()const { return node_coords_; }
        inline const auto& face_centroids()const { ASSERT(!mesh_changed_); return face_centroids_; }
        inline const auto& face_areas    ()const { ASSERT(!mesh_changed_); return face_area_; }
        inline const auto& face_normals  ()const { ASSERT(!mesh_changed_); return face_normal_; }
        inline const auto& face_types    ()const { return face_types_; }

        inline auto topo ()const { ASSERT(!mesh_changed_); return topo_; }
        inline auto shape()const { ASSERT(!mesh_changed_); return shape_; }

        inline bool is_high_order()const { ASSERT(!mesh_changed_); return is_high_order_; }

        inline auto& get_fields()const { return fields_; }

        Field&       get_field(const char* field_name);
        const Field& get_field(const char* field_name)const;

        friend class Communicator;
        friend class Interpolator;
        friend class Application;

    private:
        void register_field(const FieldInfo& fd);

    private:
        IndexSet         nodes_;
        MeshConnectivity face_nodes_, node_faces_;
        vvec             node_coords_;
        vvec             face_centroids_;
        ivec             face_types_;
        dvec             face_area_;
        vvec             face_normal_;

        //---------------------------------------------------
        // Following Data will be ignored by communicator.
        //---------------------------------------------------

        ZoneTopo         topo_ { ZT_POINTS };  //! update in add_face()
        ZoneShape        shape_{ ZS_GENERAL }; //! update in compute_metics()
        vec3             coord_min_{ 0,0,0 };
        vec3             coord_max_{ 0,0,0 };

        //!@note THe field values will be ignored by communicator.
        
        bool             mesh_changed_{ false }; //! 
        bool             is_high_order_{ false };//! Is all faces are high order element?

        Fields                   fields_;

        KDTree<double, 3, int_l> kdtree_;

        TinyMatrix<double, 4> xps_tm_;     //! used to transform point to local CS. update in compute_metics()
        dvec                  xps_coords_; //! 
        DynamicMatrix         xps_ts_inv_; //!
        bool                  xps_computed_{ false };

        std::vector<int>      ibuffer_; //! used for inverse xps_ts_inv_.
        dvec                  dbuffer_; //! used for computing XPS interpolation coefficient.
        vvec                  local_xps_points_;
    };
}
