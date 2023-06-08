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
//! @file       EasyFsi.cpp
//!             The implement of C API of EasyFsi library.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <string>

#include "KDTree.hpp"
#include "Application.hpp"
#include "SocketCommunicator.hpp"
#include "FluentCommunicator.hpp"
#include "MPICommunicator.hpp"

#include "EasyFsi.h"

//extern "C" void set_allocator(func_allocate falloc, func_deallocate fdealloc)
//{
//    EasyLib::set_user_allocator(falloc, fdealloc);
//}
//
//extern "C" void* allocate(size_t size_in_bytes)
//{
//    return EasyLib::allocate(size_in_bytes);
//}
//
//extern "C" void  deallocate(void* pointer)
//{
//    EasyLib::deallocate(pointer);
//}

//--- DynamicVector

using DV = EasyLib::DynamicVector;

extern "C" DynamicVector * dv_new(int size, double initial_value)
{
    return reinterpret_cast<DynamicVector*>(new DV(size, initial_value));
}
extern "C" void    dv_delete(DynamicVector** p_dv)
{
    if (p_dv && *p_dv) {
        delete reinterpret_cast<DV*>(*p_dv);
        *p_dv = nullptr;
    }
}
extern "C" void    dv_clear(DynamicVector * dv)
{
    if (dv)reinterpret_cast<DV*>(dv)->clear();
}
extern "C" void    dv_resize(DynamicVector * dv, int size, double initial_value)
{
    if (dv)reinterpret_cast<DV*>(dv)->resize(size, initial_value);
}
extern "C" void    dv_push_back(DynamicVector * dv, double value)
{
    if (dv)reinterpret_cast<DV*>(dv)->push_back(value);
}
extern "C" int     dv_size(const DynamicVector * dv)
{
    return dv ? (int)(reinterpret_cast<const DV*>(dv)->size()) : 0;
}
extern "C" double* dv_data(DynamicVector * dv)
{
    return dv ? reinterpret_cast<DV*>(dv)->data() : nullptr;
}
extern "C" void    dv_fill(DynamicVector * dv, double value)
{
    if (dv)reinterpret_cast<DV*>(dv)->fill(value);
}
extern "C" void    dv_copy_from(const DynamicVector* dv_src, DynamicVector* dv_des)
{
    if (dv_des)reinterpret_cast<DV*>(dv_des)->copy_elements(*reinterpret_cast<const DV*>(dv_src));
}
extern "C" void    dv_swap(DynamicVector * dvA, DynamicVector * dvB)
{
    if (dvA)reinterpret_cast<DV*>(dvA)->swap_elements(*reinterpret_cast<DV*>(dvB));
}
extern "C" double  dv_norm(const DynamicVector * dv)
{
    return dv ? reinterpret_cast<const DV*>(dv)->norm() : 0;
}
extern "C" double  dv_norm_sq(const DynamicVector* dv)
{
    return dv ? reinterpret_cast<const DV*>(dv)->norm_sq() : 0;
}
extern "C" double  dv_min(const DynamicVector* dv)
{
    return dv ? reinterpret_cast<const DV*>(dv)->min() : 0;
}
extern "C" double  dv_max(const DynamicVector* dv)
{
    return dv ? reinterpret_cast<const DV*>(dv)->max() : 0;
}
extern "C" double  dv_mean(const DynamicVector* dv)
{
    return dv ? reinterpret_cast<const DV*>(dv)->mean() : 0;
}
extern "C" double  dv_dot(const DynamicVector * dvA, const DynamicVector * dvB)
{
    return (dvA && dvB) ? dot(*reinterpret_cast<const DV*>(dvA), *reinterpret_cast<const DV*>(dvB)) : 0;
}
extern "C" double  dv_get(const DynamicVector * dv, int i)
{
    return reinterpret_cast<const DV*>(dv)->operator[](i);
}
extern "C" void    dv_set(DynamicVector * dv, int i, double value)
{
    if (dv)reinterpret_cast<DV*>(dv)->operator[](i) = value;
}

//--- DynamicMatrix

using DM = EasyLib::DynamicMatrix;

extern "C" DynamicMatrix * dm_new(int nrow, int ncol, double initial_value)
{
    return reinterpret_cast<DynamicMatrix*>(new DM(nrow, ncol, initial_value));
}
extern "C" void dm_delete(DynamicMatrix** p_dm)
{
    if (p_dm && *p_dm) {
        delete reinterpret_cast<DM*>(*p_dm);
        *p_dm = nullptr;
    }
}
extern "C" void dm_clear(DynamicMatrix * dm)
{
    if (dm)reinterpret_cast<DM*>(dm)->clear();
}
extern "C" void dm_resize(DynamicMatrix * dm, int nrow, int ncol)
{
    if (dm)reinterpret_cast<DM*>(dm)->resize(nrow, ncol);
}
extern "C" int  dm_nrow(const DynamicMatrix * dm)
{
    return dm ? (int)reinterpret_cast<const DM*>(dm)->nrow() : 0;
}
extern "C" int  dm_ncol(const DynamicMatrix * dm)
{
    return dm ? (int)reinterpret_cast<const DM*>(dm)->ncol() : 0;
}
extern "C" int  dm_numel(const DynamicMatrix * dm)
{
    return dm ? (int)reinterpret_cast<const DM*>(dm)->numel() : 0;
}
extern "C" void dm_fill(DynamicMatrix * dm, double value)
{
    if (dm)reinterpret_cast<DM*>(dm)->fill(value);
}
extern "C" void dm_copy_from(const DynamicMatrix * dm_src, DynamicMatrix * dm_des)
{
    if (dm_src && dm_des)reinterpret_cast<DM*>(dm_des)->copy_elements(*reinterpret_cast<const DM*>(dm_src));
}
//extern "C" void dm_swap(DynamicMatrix * dmA, DynamicMatrix * dmB)
//{
//    if (dmA && dmB)reinterpret_cast<DM*>(dmA)->swap();
//}
extern "C" double* dm_data(DynamicMatrix * dm)
{
    return dm ? reinterpret_cast<DM*>(dm)->data() : nullptr;
}
extern "C" double  dm_get(const DynamicMatrix * dm, int i, int j)
{
    return dm ? reinterpret_cast<const DM*>(dm)->at(i, j) : 0;
}
extern "C" void    dm_set(DynamicMatrix * dm, int i, int j, double value)
{
    if (dm)reinterpret_cast<DM*>(dm)->at(i, j) = value;
}
extern "C" void dm_identity(DynamicMatrix * dm)
{
    if (dm)reinterpret_cast<DM*>(dm)->identity();
}
extern "C" void dm_inverse(DynamicMatrix * dm, int* singular)
{
    *singular = (dm && reinterpret_cast<DM*>(dm)->inverse()) ? 1 : 0;
}
extern "C" void dm_apply_dv(const DynamicMatrix * dm, const DynamicVector * x, DynamicVector * y)
{
    if (dm && x && y)reinterpret_cast<const DM*>(dm)->apply(*reinterpret_cast<const DV*>(x), *reinterpret_cast<DV*>(y));
}
extern "C" void dm_apply_dm(const DynamicMatrix * dm, const DynamicMatrix * x, DynamicMatrix * y)
{
    if (dm && x && y)reinterpret_cast<const DM*>(dm)->apply(*reinterpret_cast<const DM*>(x), *reinterpret_cast<DM*>(y));
}
extern "C" void dm_apply_add_dv(const DynamicMatrix * dm, const DynamicVector * x, DynamicVector * y)
{
    if (dm && x && y)reinterpret_cast<const DM*>(dm)->apply_add(*reinterpret_cast<const DV*>(x), *reinterpret_cast<DV*>(y));
}
extern "C" void dm_apply_add_dm(const DynamicMatrix * dm, const DynamicMatrix * x, DynamicMatrix * y)
{
    if (dm && x && y)reinterpret_cast<const DM*>(dm)->apply_add(*reinterpret_cast<const DM*>(x), *reinterpret_cast<DM*>(y));
}

//--- IndexSet

using IS = EasyLib::IndexSet;
extern "C" IndexSet * is_new()
{
    return reinterpret_cast<IndexSet*>(new IS());
}
extern "C" void      is_delete(IndexSet** p_is)
{
    if (p_is && *p_is) {
        delete reinterpret_cast<IS*>(*p_is);
        *p_is = nullptr;
    }
}
extern "C" void      is_clear(IndexSet* is)
{
    if (is)reinterpret_cast<IS*>(is)->clear();
}
extern "C" int_l     is_add(IndexSet* is, int_g unique_id)
{
    return is ? reinterpret_cast<IS*>(is)->add(unique_id) : EasyLib::invalid_id;
}
extern "C" int       is_contains(const IndexSet* is, int_g unique_id)
{
    return is ? reinterpret_cast<const IS*>(is)->contains(unique_id) : false;
}
extern "C" int_g     is_l2g(const IndexSet* is, int_l id)
{
    return is ? reinterpret_cast<const IS*>(is)->l2g(id) : EasyLib::invalid_id;
}
extern "C" int_l     is_g2l(const IndexSet* is, int_g unique_id)
{
    return is ? reinterpret_cast<const IS*>(is)->g2l(unique_id) : EasyLib::invalid_id;
}
extern "C" int_l     is_size(const IndexSet* is)
{
    return is ? reinterpret_cast<const IS*>(is)->size() : 0;
}
extern "C" const int_g* is_glist(const IndexSet* is)
{
    return is ? reinterpret_cast<const IS*>(is)->data() : nullptr;
}

//--- KdTree

using KDT = EasyLib::KDTree<double, 3, int_l>;
extern "C" KdTree * kdt_new(const double* coords, int npts, int persistent)
{
    return reinterpret_cast<KdTree*>(new KDT(coords, npts, persistent));
}
extern "C" void    kdt_delete(KdTree** p_kdt)
{
    if (p_kdt && *p_kdt) {
        delete reinterpret_cast<KDT*>(*p_kdt);
        *p_kdt = nullptr;
    }
}
extern "C" void    kdt_clear(KdTree* kdt)
{
    if (kdt)reinterpret_cast<KDT*>(kdt)->clear();
}
extern "C" void    kdt_create(KdTree* kdt, const double* coords, int npts, int persistent)
{
    if (kdt)reinterpret_cast<KDT*>(kdt)->create(coords, npts, persistent);
}
extern "C" int     kdt_size(const KdTree* kdt)
{
    return kdt ? reinterpret_cast<const KDT*>(kdt)->size() : 0;
}
extern "C" int     kdt_search(const KdTree* kdt, const double* q, int n_query, int* ids, double* d2_sq)
{
    return kdt ? reinterpret_cast<const KDT*>(kdt)->search(q, n_query, ids, d2_sq) : 0;
}
extern "C" const double* kdt_coords(const KdTree* kdt)
{
    return kdt ? reinterpret_cast<const KDT*>(kdt)->data() : nullptr;
}

//--- MeshConnectivity

using MC = EasyLib::MeshConnectivity;

extern "C" MeshConnectivity * mc_new()
{
    return reinterpret_cast<MeshConnectivity*>(new MC());
}
extern "C" void mc_delete(MeshConnectivity * *p_mc)
{
    if (p_mc && *p_mc) {
        delete reinterpret_cast<MC*>(*p_mc);
        *p_mc = nullptr;
    }
}
extern "C" void mc_clear(MeshConnectivity * mc)
{
    if (mc)reinterpret_cast<MC*>(mc)->clear();
}
extern "C" void mc_reserve(MeshConnectivity * mc, int_l max_nrow, int_l max_ndata)
{
    if (mc)reinterpret_cast<MC*>(mc)->reserve(max_nrow, max_ndata);
}
extern "C" void mc_push_back(MeshConnectivity * mc, int n, const int_l * data)
{
    if (mc)reinterpret_cast<MC*>(mc)->push_back(n, data);
}
extern "C" const int_l * mc_ia(const MeshConnectivity * mc)
{
    return mc ? reinterpret_cast<const MC*>(mc)->ia().data() : nullptr;
}
extern "C" const int_l * mc_ja(const MeshConnectivity * mc)
{
    return mc ? reinterpret_cast<const MC*>(mc)->ja().data() : nullptr;
}
extern "C" int_l mc_nrow(const MeshConnectivity * mc)
{
    return mc ? reinterpret_cast<const MC*>(mc)->nrow() : 0;
}
extern "C" int_l mc_ndata(const MeshConnectivity * mc)
{
    return mc ? reinterpret_cast<const MC*>(mc)->ndata() : 0;
}
extern "C" int_l mc_row_size(const MeshConnectivity * mc, int_l row)
{
    return mc ? reinterpret_cast<const MC*>(mc)->ndata(row) : 0;
}
extern "C" int_l mc_row_data(const MeshConnectivity * mc, int_l row, int_l idata)
{
    return mc ? reinterpret_cast<const MC*>(mc)->operator[](row).operator[](idata) : EasyLib::invalid_id;
}


//--- Boundary

using BD = EasyLib::Boundary;

extern "C" Boundary * bd_new()
{
    return reinterpret_cast<Boundary*>(new BD());
}
extern "C" void  bd_delete(Boundary * *p_bd)
{
    if (p_bd && *p_bd) {
        delete reinterpret_cast<BD*>(*p_bd);
        *p_bd = nullptr;
    }
}
extern "C" void  bd_clear(Boundary * bd)
{
    if (bd)reinterpret_cast<BD*>(bd)->clear();
}
extern "C" void  bd_set_user_id(Boundary * bd, int id)
{
    if (bd)reinterpret_cast<BD*>(bd)->set_user_id(id);
}
extern "C" int   bd_get_user_id(const Boundary * bd)
{
    return bd ? reinterpret_cast<const BD*>(bd)->user_id() : -1;
}
extern "C" void  bd_reserve(Boundary * bd, int_l max_node, int_l max_face, int_l max_face_nodes)
{
    if (bd)reinterpret_cast<BD*>(bd)->reserve(max_node, max_face, max_face_nodes);
}
extern "C" int_l bd_add_node(Boundary * bd, double x, double y, double z, int_g unique_id)
{
    return bd ? reinterpret_cast<BD*>(bd)->add_node(x, y, z, unique_id) : EasyLib::invalid_id;
}
extern "C" int_l bd_add_face(Boundary * bd, FaceTopo type, int nnodes, const int_l * fnodes)
{
    return bd ? reinterpret_cast<BD*>(bd)->add_face((EasyLib::FaceTopo)type, nnodes, fnodes) : EasyLib::invalid_id;
}
extern "C" void  bd_set_face_centroid(Boundary * bd, int_l face, double cx, double cy, double cz)
{
    if (bd)reinterpret_cast<BD*>(bd)->set_face_cent(face, cx, cy, cz);
}
extern "C" void  bd_set_face_area(Boundary * bd, int_l face, double sx, double sy, double sz)
{
    if (bd)reinterpret_cast<BD*>(bd)->set_face_area(face, sx, sy, sz);
}
extern "C" void  bd_set_node_coords(Boundary * bd, int_l node, double x, double y, double z)
{
    if (bd)reinterpret_cast<BD*>(bd)->set_node_coords(node, x, y, z);
}
extern "C" int_l bd_face_num(const Boundary * bd)
{
    return bd ? reinterpret_cast<const BD*>(bd)->face_num() : 0;
}
extern "C" int_l bd_node_num(const Boundary * bd)
{
    return bd ? reinterpret_cast<const BD*>(bd)->node_num() : 0;
}
extern "C" void  bd_compute_metrics(Boundary * bd, double basied_angle_deg/* = 5.0*/)
{
    if (bd)reinterpret_cast<BD*>(bd)->compute_metics(basied_angle_deg);
}
extern "C" const double* bd_face_normal(const Boundary * bd, int_l face)
{
    return bd ? reinterpret_cast<const BD*>(bd)->face_normals().at(face).data() : nullptr;
}
extern "C" double        bd_face_area(const Boundary * bd, int_l face)
{
    return bd ? reinterpret_cast<const BD*>(bd)->face_areas().at(face) : 0;
}
extern "C" const double* bd_face_centroid(const Boundary * bd, int_l face)
{
    return bd ? reinterpret_cast<const BD*>(bd)->face_centroids().at(face).data() : nullptr;
}
extern "C" FaceTopo      bd_face_type(const Boundary* bd, int_l face)
{
    return bd ? (FaceTopo)reinterpret_cast<const BD*>(bd)->face_types().at(face) : BAR2;
}
extern "C" const double* bd_node_coords(const Boundary * bd, int_l node)
{
    return bd ? reinterpret_cast<const BD*>(bd)->node_coords().at(node).data() : nullptr;
}
extern "C" int_g bd_node_l2g(const Boundary* bd, int_l node)
{
    return bd ? reinterpret_cast<const BD*>(bd)->nodes().l2g(node) : -1;
}
extern "C" int_l bd_node_g2l(const Boundary * bd, int_g node)
{
    if (bd) {
        auto& x = reinterpret_cast<const BD*>(bd)->nodes();
        auto res = x.find(node);
        return res.first ? res.second : -1;
    }
    else
        return -1;
}
extern "C" const IndexSet * bd_nodes(const Boundary * bd)
{
    return bd ? (IndexSet*)&reinterpret_cast<const BD*>(bd)->nodes() : nullptr;
}
extern "C" const MeshConnectivity * bd_face_nodes(const Boundary * bd)
{
    return bd ? (MeshConnectivity*)&reinterpret_cast<const BD*>(bd)->face_nodes() : nullptr;
}
extern "C" const KdTree * bd_kdtree(const Boundary * bd)
{
    return bd ? (const KdTree*)&reinterpret_cast<const BD*>(bd)->kdtree() : nullptr;
}
extern "C" void  bd_read_gmsh(Boundary * bd, const char* file)
{
    if (bd)reinterpret_cast<BD*>(bd)->load_gmsh(file);
}

//--- Communicator

using CM = EasyLib::Communicator;
extern "C" Communicator * cm_socket_new(int as_master, int np, const char* master_ip, int master_port)
{
    char arg_np[16] = {'\0'};
    char arg_port[16] = { '\0' };

    sprintf_s(arg_np, "%d", np);
    sprintf_s(arg_port, "%d", master_port);

    const char* args[] = {
        "\0",
        as_master ? "-master" : "\0",
        "-np",
        arg_np,
        "-ip",
        master_ip ? master_ip : "\0",
        "-port",
        arg_port
    };

    auto cm = new EasyLib::SocketCommunicator();
    cm->init(8, args);
    return reinterpret_cast<Communicator*>((CM*)cm);
}
extern "C" Communicator * cm_mpi_new(int mpi_comm, int rank, int size)
{
    return reinterpret_cast<Communicator*>((CM*)new EasyLib::MPICommunicator(mpi_comm, rank, size));
}
extern "C" Communicator * cm_fluent_new(int myid, int np, func_MPT_csend * csend, func_MPT_crecv * crecv)
{
    return reinterpret_cast<Communicator*>(
        (CM*)new EasyLib::FluentCommunicator(
            myid, np,
            (EasyLib::FluentCommunicator::func_MPT_csend*)csend,
            (EasyLib::FluentCommunicator::func_MPT_crecv*)crecv
        ));
}
extern "C" void cm_set_constant(Communicator* cm, const char* name, int value)
{
    reinterpret_cast<CM*>(cm)->set_constant(name, value);
}
extern "C" void cm_set_pointer(Communicator* cm, const char* name, void* value)
{
    reinterpret_cast<CM*>(cm)->set_constant(name, value);
}
extern "C" void cm_set_function(Communicator* cm, const char* name, void* value)
{
    reinterpret_cast<CM*>(cm)->set_function(name, value);
}
extern "C" void cm_delete(Communicator * *p_cm)
{
    if (p_cm && *p_cm) {
        delete reinterpret_cast<CM*>(*p_cm);
        *p_cm = nullptr;
    }
}
extern "C" int cm_rank(Communicator * cm)
{
    return cm ? reinterpret_cast<CM*>(cm)->rank() : -1;
}
extern "C" int cm_size(Communicator * cm)
{
    return cm ? reinterpret_cast<CM*>(cm)->rank() : 01;
}
extern "C" void cm_send_int16(Communicator * cm, const int16_t * data, int count, int dest_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->send(data, count, dest_rank, tag);
}
extern "C" void cm_send_int32(Communicator * cm, const int32_t * data, int count, int dest_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->send(data, count, dest_rank, tag);
}
extern "C" void cm_send_int64(Communicator * cm, const int64_t * data, int count, int dest_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->send(data, count, dest_rank, tag);
}
extern "C" void cm_send_double(Communicator * cm, const double* data, int count, int dest_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->send(data, count, dest_rank, tag);
}
extern "C" void cm_send_float(Communicator * cm, const float* data, int count, int dest_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->send(data, count, dest_rank, tag);
}
extern "C" void cm_send_char(Communicator * cm, const char* data, int count, int dest_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->send(data, count, dest_rank, tag);
}
extern "C" void cm_recv_int16(Communicator * cm, int16_t * data, int count, int src_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->recv(data, count, src_rank, tag);
}
extern "C" void cm_recv_int32(Communicator * cm, int32_t * data, int count, int src_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->recv(data, count, src_rank, tag);
}
extern "C" void cm_recv_int64(Communicator * cm, int64_t * data, int count, int src_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->recv(data, count, src_rank, tag);
}
extern "C" void cm_recv_double(Communicator * cm, double* data, int count, int src_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->recv(data, count, src_rank, tag);
}
extern "C" void cm_recv_float(Communicator * cm, float* data, int count, int src_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->recv(data, count, src_rank, tag);
}
extern "C" void cm_recv_char(Communicator * cm, char* data, int count, int src_rank, int tag)
{
    if (cm)reinterpret_cast<CM*>(cm)->recv(data, count, src_rank, tag);
}


//--- DistributedBoundary

using DBD = EasyLib::DistributedBoundary;

extern "C" DistributedBoundary* dbd_new(Boundary * local_bd, Communicator * intra_comm, int intra_root_rank)
{
    auto db = new DBD();
    db->assemble(*reinterpret_cast<BD*>(local_bd), *reinterpret_cast<CM*>(intra_comm), intra_root_rank);
    return reinterpret_cast<DistributedBoundary*>(db);
}
extern "C" void dbd_delete(DistributedBoundary * *p_dbd)
{
    if (p_dbd && *p_dbd) {
        delete reinterpret_cast<DBD*>(*p_dbd);
        *p_dbd = nullptr;
    }
}
extern "C" void dbd_clear(DistributedBoundary * dbd)
{
    if (dbd)reinterpret_cast<DBD*>(dbd)->clear();
}
extern "C" void dbd_assemble(DistributedBoundary * dbd, Boundary * local_bd, Communicator * intra_comm, int intra_root_rank)
{
    if (dbd)reinterpret_cast<DBD*>(dbd)->assemble(
        *reinterpret_cast<BD*>(local_bd),
        *reinterpret_cast<CM*>(intra_comm),
        intra_root_rank
    );
}
extern "C" const Boundary * dbd_full_bound(const DistributedBoundary * dbd)
{
    return dbd ? (const Boundary*) &(reinterpret_cast<const DBD*>(dbd)->full_boundary()) : nullptr;
}
extern "C" const Boundary * dbd_local_bound(const DistributedBoundary * dbd)
{
    return dbd ? (const Boundary*)&(reinterpret_cast<const DBD*>(dbd)->local_boundary()) : nullptr;
}
extern "C" void dbd_gather_node_fields(DistributedBoundary * dbd, int nfields, const double* local_fields, double* global_fields)
{
    if (dbd)reinterpret_cast<DBD*>(dbd)->gather_node_fields(nfields, local_fields, global_fields);
}
extern "C" void dbd_gather_face_fields(DistributedBoundary * dbd, int nfields, const double* local_fields, double* global_fields)
{
    if (dbd)reinterpret_cast<DBD*>(dbd)->gather_face_fields(nfields, local_fields, global_fields);
}
extern "C" void dbd_scatter_node_fields(DistributedBoundary * dbd, int nfields, const double* global_fields, double* local_fields)
{
    if (dbd)reinterpret_cast<DBD*>(dbd)->scatter_node_fields(nfields, global_fields, local_fields);
}
extern "C" void dbd_scatter_face_fields(DistributedBoundary * dbd, int nfields, const double* global_fields, double* local_fields)
{
    if (dbd)reinterpret_cast<DBD*>(dbd)->scatter_face_fields(nfields, global_fields, local_fields);
}
extern "C" void dbd_accumulate_node_fields(DistributedBoundary * dbd, int nfields, const double* local_fields, double* global_fields)
{
    if (dbd)reinterpret_cast<DBD*>(dbd)->accumulate_node_fields(nfields, local_fields, global_fields);
}
extern "C" void dbd_accumulate_face_fields(DistributedBoundary * dbd, int nfields, const double* local_fields, double* global_fields)
{
    if (dbd)reinterpret_cast<DBD*>(dbd)->accumulate_face_fields(nfields, local_fields, global_fields);
}

//--- Interpolator

using IT = EasyLib::Interpolator;
extern "C" Interpolator* it_new()
{
    return reinterpret_cast<Interpolator*>(new IT);
}
extern "C" void it_delete(Interpolator** it)
{
    if (it && *it) {
        delete reinterpret_cast<IT*>(*it);
        *it = nullptr;
    }
}
extern "C" void it_clear(Interpolator * it)
{
    if (it)reinterpret_cast<IT*>(it)->clear();
}
extern "C" int it_add_source_boundary(Interpolator* it, Boundary* bd)
{
    return (it && bd) ? reinterpret_cast<IT*>(it)->add_source_boundary(*reinterpret_cast<BD*>(bd)) : -1;
}
extern "C" int it_add_target_boundary(Interpolator* it, Boundary* bd)
{
    return (it && bd) ? reinterpret_cast<IT*>(it)->add_target_boundary(*reinterpret_cast<BD*>(bd)) : -1;
}
extern "C" void it_compute_interp_coeff(Interpolator* it, InterpMethod method, int max_donor)
{
    if (it)reinterpret_cast<IT*>(it)->compute_interp_coeff((EasyLib::InterpolationMethod)method, max_donor);
}
extern "C" void it_save_coefficients(Interpolator* it, const char* file)
{
    if (it)reinterpret_cast<IT*>(it)->save_coefficients(file);
}
extern "C" void it_load_coefficients(Interpolator* it, const char* file)
{
    if (it)reinterpret_cast<IT*>(it)->load_coefficients(file);
}
extern "C" void it_interp_all_dofs_s2t(Interpolator* it)
{
    if (it)reinterpret_cast<IT*>(it)->interp_all_dofs_s2t();
}
extern "C" void it_interp_all_load_t2s(Interpolator* it)
{
    if (it)reinterpret_cast<IT*>(it)->interp_all_load_t2s();
}
extern "C" void it_interp_node_dofs_s2t(Interpolator * it, int ndof, const double** src_node_dofs, double** des_node_dofs)
{
    if (it)reinterpret_cast<IT*>(it)->interp_node_dofs_s2t(ndof, src_node_dofs, des_node_dofs);
}
extern "C" void it_interp_face_dofs_s2t(Interpolator * it, int ndof, const double** src_node_dofs, double** des_face_dofs)
{
    if (it)reinterpret_cast<IT*>(it)->interp_face_dofs_s2t(ndof, src_node_dofs, des_face_dofs);
}
extern "C" void it_interp_node_loads_t2s(Interpolator * it, int nload, double** src_node_load, const double** des_node_load, bool fill_src_zeros_first/* = true*/)
{
    if (it)reinterpret_cast<IT*>(it)->interp_node_load_t2s(nload, src_node_load, des_node_load, fill_src_zeros_first);
}
extern "C" void it_interp_face_loads_t2s(Interpolator * it, int nload, double** src_node_load, const double** des_face_load, bool fill_src_zeros_first/* = true*/)
{
    if (it)reinterpret_cast<IT*>(it)->interp_face_load_t2s(nload, src_node_load, des_face_load, fill_src_zeros_first);
}

//--- Application

using APP = EasyLib::Application;

extern "C" Application * app_new(const char* name, Communicator * intra_comm, int root)
{
    return reinterpret_cast<Application*>(new APP(name, *reinterpret_cast<CM*>(intra_comm), root));
}
extern "C" void app_delete(Application** p_app)
{
    if (p_app && *p_app) {
        delete reinterpret_cast<APP*>(*p_app);
        *p_app = nullptr;
    }
}
extern "C" void app_clear(Application * app)
{
    if (app)reinterpret_cast<APP*>(app)->clear();
}
extern "C" Boundary * app_add_boundary(Application * app)
{
    return app ? reinterpret_cast<Boundary*>(&(reinterpret_cast<APP*>(app)->add_coupled_boundary())) : nullptr;
}
extern "C" int app_boundary_num(const Application* app)
{
    return app ? reinterpret_cast<const APP*>(app)->boundary_num() : 0;
}
extern "C" Boundary* app_get_boundary(Application* app, int ib)
{
    return app ? reinterpret_cast<Boundary*>(reinterpret_cast<APP*>(app)->boundary(ib)) : nullptr;
}
extern "C" void app_register_field(Application * app, const char* name, int ncomp, FieldLocation location, FieldIO iotype, const char* units)
{
    if (app)reinterpret_cast<APP*>(app)->register_field(
        name, ncomp,
        (EasyLib::FieldLocation)location,
        (EasyLib::FieldIO)iotype,
        units);
}
extern "C" void app_set_field_func(Application* app, get_boundary_field_function getter, set_boundary_field_function setter)
{
    if (app) {
        reinterpret_cast<APP*>(app)->set_field_function(
            (EasyLib::get_boundary_field_function)getter,
            (EasyLib::set_boundary_field_function)setter
        );
    }
}
extern "C" void app_start_coupling(Application * app, Communicator * inter_comm)
{
    if (app)reinterpret_cast<APP*>(app)->start_coupling(*reinterpret_cast<CM*>(inter_comm));
}
extern "C" void app_exchange_solu(Application* app, double time, void* user_data)
{
    if (app)reinterpret_cast<APP*>(app)->exchange_solution(time, user_data);
}
extern "C" void app_stop_coupling(Application * app)
{
    if (app)reinterpret_cast<APP*>(app)->stop_coupling();
}
