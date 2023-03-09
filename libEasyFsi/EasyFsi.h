#pragma once
#ifndef _EASYFSI_H_
#define _EASYFSI_H_
#include <stdint.h>

typedef int       int_l;
typedef long long int_g;

typedef struct DynamicVector       DynamicVector;
typedef struct DynamicMatrix       DynamicMatrix;
typedef struct IndexSet            IndexSet;
typedef struct KdTree              KdTree;
typedef struct MeshConnectivity    MeshConnectivity;
typedef struct Boundary            Boundary;
typedef struct DistributedBoundary DistributedBoundary;
typedef struct Communicator        Communicator;
typedef struct Application         Application;

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

typedef void(__stdcall* get_boundary_field_function)(const Boundary* bd, const char* name, int ncomp, FieldLocation loc,       double* data, void* user_data);
typedef void(__stdcall* set_boundary_field_function)(const Boundary* bd, const char* name, int ncomp, FieldLocation loc, const double* data, void* user_data);

typedef void(*func_MPT_csend)(int mpid, void* data, unsigned int n, int data_type, int tag, const char* file, int line);
typedef int (*func_MPT_crecv)(int mpid, void* data, unsigned int n, int data_type, int tag, const char* file, int line);

#ifdef __cplusplus
extern "C" {
#endif

    //--- DynamicVector

    DynamicVector* dv_new(int size, double initial_value);
    void    dv_delete   (DynamicVector** p_dv);
    void    dv_clear    (DynamicVector* dv);
    void    dv_resize   (DynamicVector* dv, int size, double initial_value);
    void    dv_push_back(DynamicVector* dv, double value);
    int     dv_size     (const DynamicVector* dv);
    double* dv_data     (DynamicVector* dv);
    void    dv_fill     (DynamicVector* dv, double value);
    void    dv_copy_from(const DynamicVector* dv_src, DynamicVector* dv_des);
    void    dv_swap     (DynamicVector* dvA, DynamicVector* dvB);
    double  dv_norm     (const DynamicVector* dv);
    double  dv_norm_sq  (const DynamicVector* dv);
    double  dv_min      (const DynamicVector* dv);
    double  dv_max      (const DynamicVector* dv);
    double  dv_mean     (const DynamicVector* dv);
    double  dv_dot      (const DynamicVector* dvA, const DynamicVector* dvB);
    double  dv_get      (const DynamicVector* dv, int i);
    void    dv_set      (DynamicVector* dv, int i, double value);

    //--- DynamicMatrix

    DynamicMatrix* dm_new(int nrow, int ncol, double initial_value);
    void dm_delete(DynamicMatrix** p_dm);
    void dm_clear (DynamicMatrix* dm);
    void dm_resize(DynamicMatrix* dm, int nrow, int ncol);
    int  dm_nrow(const DynamicMatrix* dm);
    int  dm_ncol(const DynamicMatrix* dm);
    int  dm_numel(const DynamicMatrix* dm);
    void dm_fill(DynamicMatrix* dm, double value);
    void dm_copy_from(const DynamicMatrix* dm_src, DynamicMatrix* dm_des);
    double* dm_data(DynamicMatrix* dm);
    double  dm_get(const DynamicMatrix* dm, int i, int j);
    void    dm_set(DynamicMatrix* dm, int i, int j, double value);
    void dm_identity(DynamicMatrix* dm);
    void dm_inverse(DynamicMatrix* dm, int* singular);
    void dm_apply_dv(const DynamicMatrix* dm, const DynamicVector* x, DynamicVector* y);
    void dm_apply_dm(const DynamicMatrix* dm, const DynamicMatrix* x, DynamicMatrix* y);
    void dm_apply_add_dv(const DynamicMatrix* dm, const DynamicVector* x, DynamicVector* y);
    void dm_apply_add_dm(const DynamicMatrix* dm, const DynamicMatrix* x, DynamicMatrix* y);

    //--- IndexSet

    IndexSet* is_new();
    void      is_delete(IndexSet** p_is);
    void      is_clear(IndexSet* is);
    int_l     is_add(IndexSet* is, int_g unique_id);
    int       is_contains(const IndexSet* is, int_g unique_id);
    int_g     is_l2g(const IndexSet* is, int_l id);
    int_l     is_g2l(const IndexSet* is, int_g unique_id);
    int_l     is_size(const IndexSet* is);
    const int_g* is_glist(const IndexSet* is);

    //--- KdTree

    KdTree* kdt_new(const double* coords, int npts, int persistent);
    void    kdt_delete(KdTree** p_kdt);
    void    kdt_clear(KdTree* kdt);
    void    kdt_create(KdTree* kdt, const double* coords, int npts, int persistent);
    int     kdt_size(const KdTree* kdt);
    int     kdt_search(const KdTree* kdt, const double* q, int n_query, int* ids, double* d2_sq);
    const double* kdt_coords(const KdTree* kdt);

    //--- MeshConnectivity

    MeshConnectivity* mc_new();
    void mc_delete(MeshConnectivity** p_mc);
    void mc_clear(MeshConnectivity* mc);
    void mc_reserve(MeshConnectivity* mc, int_l max_nrow, int_l max_ndata);
    void mc_push_back(MeshConnectivity* mc, int n, const int_l* data);
    const int_l* mc_ia(const MeshConnectivity* mc);
    const int_l* mc_ja(const MeshConnectivity* mc);
    int_l mc_nrow (const MeshConnectivity* mc);
    int_l mc_ndata(const MeshConnectivity* mc);
    int_l mc_row_size(const MeshConnectivity* mc, int_l row);
    int_l mc_row_data(const MeshConnectivity* mc, int_l row, int_l idata);

    //--- Boundary

    Boundary* bd_new();
    void  bd_delete(Boundary** p_bd);
    void  bd_clear(Boundary* bd);
    void  bd_set_user_id(Boundary* bd, int id);
    int   bd_get_user_id(Boundary* bd);
    void  bd_reserve(Boundary* bd, int_l max_node, int_l max_face, int_l max_face_nodes);
    int_l bd_add_node(Boundary* bd, double x, double y, double z, int_g unique_id);
    int_l bd_add_face(Boundary* bd, FaceTopo type, int nnodes, const int_l* fnodes);
    void  bd_set_face_centroid(Boundary* bd, int_l face, double cx, double cy, double cz);
    void  bd_set_face_area    (Boundary* bd, int_l face, double sx, double sy, double sz);
    void  bd_set_node_coords  (Boundary* bd, int_l node, double x, double y, double z);
    int_l bd_face_num(const Boundary* bd);
    int_l bd_node_num(const Boundary* bd);
    void  bd_compute_metrics(Boundary* bd, double basied_angle_deg = 5.0);
    const double* bd_face_normal  (const Boundary* bd, int_l face);
    double        bd_face_area    (const Boundary* bd, int_l face);
    const double* bd_face_centroid(const Boundary* bd, int_l face);
    const double* bd_node_coords  (const Boundary* bd, int_l node);

    const IndexSet* bd_nodes(const Boundary* bd);
    const MeshConnectivity* bd_face_nodes(const Boundary* bd);
    const KdTree* bd_kdtree(const Boundary* bd);
    void  bd_read_gmsh(Boundary* bd, const char* file);

    //--- Communicator

    Communicator* cm_socket_new(int as_master, int np, const char* master_ip, int master_port);
    Communicator* cm_mpi_new(int mpi_comm);
    Communicator* cm_fluent_new(int myid, int np, func_MPT_csend* csend, func_MPT_crecv* crecv);
    void cm_set_constant(Communicator* cm, const char* name, int value);
    void cm_set_pointer (Communicator* cm, const char* name, void* value);
    void cm_set_function(Communicator* cm, const char* name, void* value);
    void cm_delete(Communicator** p_cm);
    int  cm_rank(Communicator* cm);
    int  cm_size(Communicator* cm);
    void cm_send_int16 (Communicator* cm, const int16_t* data, int count, int dest_rank, int tag);
    void cm_send_int32 (Communicator* cm, const int32_t* data, int count, int dest_rank, int tag);
    void cm_send_int64 (Communicator* cm, const int64_t* data, int count, int dest_rank, int tag);
    void cm_send_double(Communicator* cm, const double* data, int count, int dest_rank, int tag);
    void cm_send_float (Communicator* cm, const float* data, int count, int dest_rank, int tag);
    void cm_send_char  (Communicator* cm, const char* data, int count, int dest_rank, int tag);
    void cm_recv_int16 (Communicator* cm, int16_t* data, int count, int src_rank, int tag);
    void cm_recv_int32 (Communicator* cm, int32_t* data, int count, int src_rank, int tag);
    void cm_recv_int64 (Communicator* cm, int64_t* data, int count, int src_rank, int tag);
    void cm_recv_double(Communicator* cm, double* data, int count, int src_rank, int tag);
    void cm_recv_float (Communicator* cm, float* data, int count, int src_rank, int tag);
    void cm_recv_char  (Communicator* cm, char* data, int count, int src_rank, int tag);

    //--- DistributedBoundary

    DistributedBoundary* dbd_new(Boundary* local_bd, Communicator* intra_comm, int intra_root_rank);
    void dbd_delete(DistributedBoundary** p_dbd);
    void dbd_clear(DistributedBoundary* dbd);
    void dbd_assemble(DistributedBoundary* dbd, Boundary* local_bd, Communicator* intra_comm, int intra_root_rank);
    const Boundary* dbd_full_bound(const DistributedBoundary* dbd);
    const Boundary* dbd_local_bound(const DistributedBoundary* dbd);
    void dbd_gather_node_fields(DistributedBoundary* dbd, int nfields, const double* local_fields, double* global_fields);
    void dbd_gather_face_fields(DistributedBoundary* dbd, int nfields, const double* local_fields, double* global_fields);
    void dbd_scatter_node_fields(DistributedBoundary* dbd, int nfields, const double* global_fields, double* local_fields);
    void dbd_scatter_face_fields(DistributedBoundary* dbd, int nfields, const double* global_fields, double* local_fields);
    void dbd_accumulate_node_fields(DistributedBoundary* dbd, int nfields, const double* local_fields, double* global_fields);
    void dbd_accumulate_face_fields(DistributedBoundary* dbd, int nfields, const double* local_fields, double* global_fields);

    
    //--- Application

    Application* app_new(const char* name, Communicator* intra_comm, int root);
    void app_delete(Application** p_app);
    void app_clear(Application* app);
    void app_create(Application* app, const char* name, Communicator* intra_comm, int root);
    Boundary* app_add_boundary(Application* app);
    void app_register_field(Application* app, const char* name, int ncomp, FieldLocation location, FieldIO iotype, const char* units);
    void app_start_coupling(Application* app, Communicator* inter_comm);
    void app_exchange_solu(Application* app, get_boundary_field_function getter, set_boundary_field_function setter, double time, void* user_data);
    void app_stop_coupling(Application* app);

#ifdef __cplusplus
}
#endif

#endif /*_EASYFSI_H_*/
