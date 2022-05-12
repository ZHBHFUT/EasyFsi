#pragma once
#include <vector>
#include <set>

#include "Index.hpp"
#include "TinyVector.hpp"
#include "TinyMatrix.hpp"
#include "IndexSet.hpp"
#include "MeshConnectivity.hpp"
#include "KDTree.hpp"
#include "DynamicMatrix.hpp"

namespace EasyLib {

    //! @brief 边界网格单元形状
    typedef enum FaceTopo
    {
        FT_BAR2 = 0,//! 2节点线单元（线性）
        FT_BAR3,    //! 3节点线单元（二次）
        FT_TRI3,    //! 3节点三角形单元（线性）
        FT_TRI6,    //! 6节点三角形单元（二次）
        FT_QUAD4,   //! 4节点四边形单元（线性）
        FT_QUAD8,   //! 8节点四边形单元（二次）
        FT_POLYGON  //! 多边形单元（节点数 > 4)
    }FaceTopo;
    //! @brief 面周围节点个数
    extern const int npf[FT_POLYGON + 1];
    //! @brief 最大面周围节点个数
    inline constexpr int npf_max = 8;
    //! @brief 面阶次，1=线性单元，2=二次单元
    extern const int face_order[FT_POLYGON + 1];

    //! @brief 边界类型
    typedef enum ZoneType
    {
        ZT_POINTS = 0, //! 离散点集合
        ZT_CURVE,      //! 曲线
        ZT_SURFACE,    //! 曲面
        ZT_VOLUME      //! 体
    }ZoneType;

    //! @brief 边界形状
    typedef enum ZoneTopo
    {
        TP_POINT = 0, //! 所有点重合
        TP_COLINEAR,  //! 所有点共线
        TP_COPLANER,  //! 所有点共面
        TP_GENERAL    //! 一般三维空间分布
    }ZoneTopo;

    //! @brief 边界数据
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

        void clear();

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
        //! @param id 
        //! @param x 
        //! @param y 
        //! @param z 
        void set_node_coords(int_l id, double x, double y, double z);

        void set_face_cent(int_l face, double cx, double cy, double cz);

        void set_face_area(int_l face, double sx, double sy, double sz);

        //! @brief Compute area, normal, centroid of faces, and topo type.
        //! @note This member function should be called after add all nodes and faces, or after modifying all node coordinates.
        void compute_metics(double biased_angle_deg = 5);

        inline auto& kdtree()const { ASSERT(!mesh_changed_); return kdtree_; }

        void compute_global_xps_matrix();

        void compute_global_xps_interp_coeff(const vec3& p, std::span<double> coeff);

        void compute_local_xps_interp_coeff(const vec3& p, int max_neigh, std::span<int_l> ids, std::span<double> coeff, int& count);

        void compute_project_interp_coeff(const vec3& p, int_l(&ids)[npf_max], double(&coeff)[npf_max], int& count);

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

        inline ZoneType type()const { ASSERT(!mesh_changed_); return type_; }
        inline ZoneTopo topo()const { ASSERT(!mesh_changed_); return topo_; }

        double* allocate_node_fields(int nfields);
        double* allocate_face_fields(int nfields);
        void    delete_fields       (double** fields);

        friend class Communicator;
        friend class Interpolator;
        
    private:
        IndexSet         nodes_;
        MeshConnectivity face_nodes_, node_faces_;
        vvec             node_coords_;
        vvec             face_centroids_;
        ivec             face_types_;
        dvec             face_area_;
        vvec             face_normal_;

        ZoneType         type_{ ZT_POINTS };  //! update in add_face()
        ZoneTopo         topo_{ TP_GENERAL }; //! update in compute_metics()
        vec3             coord_min_{ 0,0,0 };
        vec3             coord_max_{ 0,0,0 };

        bool             mesh_changed_{ false };

        //---------------------------------------------------
        // data bellow will not be processed by communicator
        //

        std::set<double*>        fields_;
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
