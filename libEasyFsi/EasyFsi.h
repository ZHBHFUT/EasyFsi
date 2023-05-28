#pragma once
#ifndef _EASYFSI_H_
#define _EASYFSI_H_
#include <stdint.h>

typedef int       int_l; //! integer type used to represent local grid item id.
typedef long long int_g; //! integer type used to represent global grid item id.

typedef struct DynamicVector       DynamicVector;
typedef struct DynamicMatrix       DynamicMatrix;
typedef struct IndexSet            IndexSet;
typedef struct KdTree              KdTree;
typedef struct MeshConnectivity    MeshConnectivity;
typedef struct Boundary            Boundary;
typedef struct DistributedBoundary DistributedBoundary;
typedef struct Communicator        Communicator;
typedef struct Application         Application;

//! @brief The location of field.
enum FieldLocation
{
    NodeCentered = 0, //! data is stored on nodes 
    FaceCentered = 1, //! data is stored on faces
    CellCentered = 2  //! data is stored on cells
};

//! @brief The input/output type of field.
enum FieldIO
{
    IncomingDofs  = 0, //! field data will be treated as DOFs and received from other application.
    IncomingLoads = 1, //! field data will be treated as LOADs and received from other application.
    OutgoingDofs  = 2, //! field data will be treated as DOFs and sending to other application.
    OutgoingLoads = 3  //! field data will be treated as LOADs and sending to other application.
};

//! @brief Topology type of face on coupled boundary.
typedef enum ElementShape
{
    BAR2 = 0,//! 2-nodes linear element
    BAR3,    //! 3-nodes quadratic element
    TRI3,    //! 3-nodes linear triangle element
    TRI6,    //! 6-nodes quadratic triangle element
    QUAD4,   //! 4-node linear quadrilateral element
    QUAD8,   //! 8-node quadratic quadrilateral element
    POLYGON  //! general polygon element (node number > 4)
}FaceTopo;

//! @brief function prototype used to read boundary field data from current application.
typedef void(__stdcall* get_boundary_field_function)(const Boundary* bd, const char* name, int ncomp, FieldLocation loc,       double* data, void* user_data);
//! @brief function prototype used to write boundary field data to current application.
typedef void(__stdcall* set_boundary_field_function)(const Boundary* bd, const char* name, int ncomp, FieldLocation loc, const double* data, void* user_data);

//! @brief function prototype used to send data to other process of this application.
typedef void(__stdcall *func_MPT_csend)(int mpid, void* data, unsigned int n, int data_type, int tag, const char* file, int line);
//! @brief function prototype used to receive data from other process of this application.
typedef int (__stdcall *func_MPT_crecv)(int mpid, void* data, unsigned int n, int data_type, int tag, const char* file, int line);

//! @brief function prototype used to send data to other process of this application.
typedef int(__stdcall* func_MPI_Send)(const void* buffer, int count, int datatype, int dest,   int tag, int comm);
//! @brief function prototype used to receive data from other process of this application.
typedef int(__stdcall* func_MPI_Recv)(      void* buffer, int count, int datatype, int source, int tag, int comm);

typedef void* (*func_allocate  )(size_t nbytes);
typedef void  (*func_deallocate)(void* pointer);

#ifdef __cplusplus
extern "C" {
#endif

    void set_allocator(func_allocate falloc, func_deallocate fdealloc);

    void* allocate(size_t size_in_bytes);

    void  deallocate(void* pointer);

    //--------------------------------------------------------
    // interface of DynamicVector
    //--------------------------------------------------------
    
    //! @brief Create a new dynamic vector.
    //! @param size          The initial size of vector
    //! @param initial_value The initial value of vector
    //! @return Return a pointer of dynamic vector object.
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

    //--------------------------------------------------------
    // interface of DynamicMatrix
    //--------------------------------------------------------

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

    //--------------------------------------------------------
    // interface of IndexSet
    //--------------------------------------------------------

    IndexSet* is_new();
    void      is_delete(IndexSet** p_is);
    void      is_clear(IndexSet* is);
    int_l     is_add(IndexSet* is, int_g unique_id);
    int       is_contains(const IndexSet* is, int_g unique_id);
    int_g     is_l2g(const IndexSet* is, int_l id);
    int_l     is_g2l(const IndexSet* is, int_g unique_id);
    int_l     is_size(const IndexSet* is);
    const int_g* is_glist(const IndexSet* is);

    //--------------------------------------------------------
    // interface of KdTree
    //--------------------------------------------------------

    KdTree* kdt_new(const double* coords, int npts, int persistent);
    void    kdt_delete(KdTree** p_kdt);
    void    kdt_clear(KdTree* kdt);
    void    kdt_create(KdTree* kdt, const double* coords, int npts, int persistent);
    int     kdt_size(const KdTree* kdt);
    int     kdt_search(const KdTree* kdt, const double* q, int n_query, int* ids, double* d2_sq);
    const double* kdt_coords(const KdTree* kdt);

    //--------------------------------------------------------
    // interface of MeshConnectivity
    //--------------------------------------------------------

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
    
    //--------------------------------------------------------
    // interface of Boundary
    //--------------------------------------------------------

    Boundary* bd_new();
    void  bd_delete(Boundary** p_bd);
    void  bd_clear(Boundary* bd);
    void  bd_set_user_id(Boundary* bd, int id);
    int   bd_get_user_id(const Boundary* bd);
    void  bd_reserve(Boundary* bd, int_l max_node, int_l max_face, int_l max_face_nodes);
    int_l bd_add_node(Boundary* bd, double x, double y, double z, int_g unique_id);
    int_l bd_add_face(Boundary* bd, ElementShape type, int nnodes, const int_l* fnodes);
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

    int_g bd_node_l2g(const Boundary* bd, int_l node);
    int_l bd_node_g2l(const Boundary* bd, int_g node);

    const IndexSet* bd_nodes(const Boundary* bd);
    const MeshConnectivity* bd_face_nodes(const Boundary* bd);
    const KdTree* bd_kdtree(const Boundary* bd);
    void  bd_read_gmsh(Boundary* bd, const char* file);

    //--------------------------------------------------------
    // interface of Communicator
    //--------------------------------------------------------

    Communicator* cm_socket_new(int as_master, int np, const char* master_ip, int master_port);
    Communicator* cm_mpi_new(int mpi_comm, int rank, int size);
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

    //--------------------------------------------------------
    // interface of DistributedBoundary
    //--------------------------------------------------------

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

    //--------------------------------------------------------
    // interfaces of Application
    //--------------------------------------------------------

    //! @brief Create a new application object.
    //! @param name        The name of this application
    //! @param intra_comm  The intra-communicator, usually is MPICommunicator, null=not exists
    //! @param root        The root process used to run FSI, ignored if \intra_comm is null.
    //! @return  Return a pointer of application object.
    Application* app_new(const char* name, Communicator* intra_comm, int root);

    //! @brief Delete an application.
    //! @param p_app  A pointer of object pointer created by \app_new.
    void app_delete(Application** p_app);

    //! @brief Remove all boundaries, communicators and other data in the application.
    //! @param app The object pointer of application created by \app_new.
    void app_clear(Application* app);

    //! @brief Add a new boundary to the application.
    //! @param app  The object pointer of application created by \app_new.
    //! @return  Return a pointer of boundary object.
    Boundary* app_add_boundary(Application* app);

    int app_boundary_num(const Application* app);

    Boundary* app_get_boundary(Application* app, int ib);

    //! @brief Register a field to the application.
    //! @param app       The object pointer of application created by \app_new.
    //! @param name      The name of the field.
    //! @param ncomp     The number of components of the field, should be positive.
    //! @param location  The data location of this field, see \FieldLocation.
    //! @param iotype    The input/output type of this field, see \FieldIO.
    //! @param units     The units of this field.
    void app_register_field(Application* app, const char* name, int ncomp, FieldLocation location, FieldIO iotype, const char* units);

    //! @brief Start coupling between applications.
    //! @param app        The object pointer of application created by \app_new.
    //! @param inter_comm The inter-communicator, usually is socket communicator.
    //! @note This call will be blocking until finishing operations bellow:
    //!    1) assemble boundary
    //!    2) synchronize all fields information between applications
    //!    3) allocate fields
    //!    4) compute interpolation coefficients.
    void app_start_coupling(Application* app, Communicator* inter_comm);

    //! @brief Exchange fields between applications.
    //! @param app       The object pointer of application created by \app_new.
    //! @param getter    A function pointer used to read outgoing fields from this application.
    //! @param setter    A function pointer used to writ incoming fields to this application.
    //! @param time      Current physical time of this application.
    //! @param user_data User data will be transfer to \getter and \setter functions, can be null.
    void app_exchange_solu(Application* app, get_boundary_field_function getter, set_boundary_field_function setter, double time, void* user_data);

    //! @brief Stop coupling of all applications and disconnect from inter-communicator.
    //! @param app  The object pointer of application created by \app_new.
    void app_stop_coupling(Application* app);

#ifdef __cplusplus
}
#endif

#endif /*_EASYFSI_H_*/
