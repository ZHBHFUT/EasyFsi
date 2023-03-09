module easyfsi
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#ifdef HAVE_EASYFSI

use iso_c_binding,only: c_int,c_float,c_double,c_long_long,c_char,c_ptr,c_funptr

implicit none

enum, bind(c)
    enumerator :: NodeCentered = 0
    enumerator :: FaceCentered = 1
    enumerator :: CellCentered = 2
end enum

enum, bind(c)
    enumerator :: IncomingDofs  = 0
    enumerator :: IncomingLoads = 1
    enumerator :: OutgoingDofs  = 2
    enumerator :: OutgoingLoads = 4
end enum

enum, bind(c)
    enumerator :: FT_BAR2 = 0
    enumerator :: FT_BAR3 = 1
    enumerator :: FT_TRI3 = 2
    enumerator :: FT_TRI6 = 3
    enumerator :: FT_QUAD4 = 4
    enumerator :: FT_QUAD8 = 5
    enumerator :: FT_POLYGON = 6
end enum

!-----------------------------------------------------
! function prototype of get_boundary_field
!-----------------------------------------------------
!
! subroutine get_boundary_field(bd,name,ncomp,location,data,user_data)
! !DEC$ ATTRIBUTES STDCALL
!     use iso_c_binding, only: c_ptr,c_char
!     type(c_ptr), intent(in), value :: bd
!     character(kind=c_char), dimension(*), intent(in) :: name
!     integer(kind=c_int), intent(in), value :: ncomp
!     integer(kind(NodeCentered)), intent(in), value :: location
!     real(kind=c_double), dimension(*), intent(out) :: data
!     type(c_ptr), intent(in), value :: user_data
! end subroutine
!

!-----------------------------------------------------
! function prototype of set_boundary_field
!-----------------------------------------------------
!
! subroutine set_boundary_field(bd,name,ncomp,location,data,user_data)
! !DEC$ ATTRIBUTES STDCALL
!     use iso_c_binding, only: c_ptr,c_char
!     type(c_ptr), intent(in), value :: bd
!     character(kind=c_char), dimension(*), intent(in) :: name
!     integer(kind=c_int), intent(in), value :: ncomp
!     integer(kind(NodeCentered)), intent(in), value :: location
!     real(kind=c_double), dimension(*), intent(in) :: data
!     type(c_ptr), intent(in), value :: user_data
! end subroutine
!

interface

function mpi_initialized_t(flag) bind(c)
    import :: c_int
    integer(c_int), intent(out) :: flag
    integer(c_int) :: mpi_initialized_t
end function mpi_initialized_t
function mpi_init_t(argc,argv) bind(c)
    import :: c_int,c_ptr
    integer(c_int), intent(in) :: argc
    type(c_ptr), intent(in), value :: argv
    integer(c_int) :: mpi_init_t
end function mpi_init_t
function mpi_comm_disconnect_t(comm) bind(c)
    import :: c_int
    integer(c_int), intent(inout) :: comm
    integer(c_int) :: mpi_comm_disconnect_t
end function mpi_comm_disconnect_t
function mpi_comm_rank_t(comm,rank) bind(c)
    import :: c_int
    integer(c_int), intent(in), value :: comm
    integer(c_int), intent(out) :: rank
    integer(c_int) :: mpi_comm_rank_t
end function mpi_comm_rank_t
function mpi_comm_size_t(comm,size) bind(c)
    import :: c_int
    integer(c_int), intent(in), value :: comm
    integer(c_int), intent(out) :: size
    integer(c_int) :: mpi_comm_size_t
end function mpi_comm_size_t
function mpi_send_t(buffer,count,datatype,dest,tag,comm) bind(c)
    import :: c_int,c_ptr
    type(c_ptr), intent(in), value :: buffer
    integer(c_int), intent(in), value :: count
    integer(c_int), intent(in), value :: datatype
    integer(c_int), intent(in), value :: dest
    integer(c_int), intent(in), value :: tag
    integer(c_int), intent(in), value :: comm
    integer(c_int) :: mpi_send_t
end function mpi_send_t
function mpi_recv_t(buffer,count,datatype,source,tag,comm,status) bind(c)
    import :: c_int,c_ptr
    type(c_ptr), intent(in), value :: buffer
    integer(c_int), intent(in), value :: count
    integer(c_int), intent(in), value :: datatype
    integer(c_int), intent(in), value :: source
    integer(c_int), intent(in), value :: tag
    integer(c_int), intent(in), value :: comm
    type(c_ptr), intent(in), value :: status
    integer(c_int) :: mpi_recv_t
end function mpi_recv_t
subroutine get_boundary_field_t(bd,fieldname,ncomp,location,data,user_data) bind(c)
    import :: c_ptr,c_int,c_char,c_double
    type(c_ptr),                   intent(in), value :: bd
    character(kind=c_char), dimension(*), intent(in) :: fieldname
    integer(kind=c_int),           intent(in), value :: ncomp
    integer(kind(NodeCentered)),   intent(in), value :: location
    real(kind=c_double),   dimension(*), intent(out) :: data
    type(c_ptr),                   intent(in), value :: user_data
end subroutine get_boundary_field_t
subroutine set_boundary_field_t(bd,fieldname,ncomp,location,fielddata,user_data) bind(c)
    import :: c_ptr,c_char,c_double,c_int
    type(c_ptr),                   intent(in), value :: bd
    character(c_char),      dimension(*), intent(in) :: fieldname
    integer(c_int),                intent(in), value :: ncomp
    integer(kind(NodeCentered)),   intent(in), value :: location
    real(kind=c_double),    dimension(*), intent(in) :: fielddata
    type(c_ptr),                   intent(in), value :: user_data
end subroutine set_boundary_field_t

!-----------------------------------------
! interfaces for application
!-----------------------------------------

function app_new(appname,intra_comm,root) bind(c,name="app_new")
    import :: c_ptr,c_int,c_char
    character(kind=c_char), dimension(*), intent(in) :: appname
    type(c_ptr), intent(in), value :: intra_comm
    integer(c_int), intent(in), value :: root
    type(c_ptr) :: app_new
end function app_new

subroutine app_delete(app) bind(c,name="app_delete")
    import :: c_ptr
    type(c_ptr), intent(inout) :: app
end subroutine app_delete

subroutine app_clear(app) bind(c,name="app_clear")
    import :: c_ptr
    type(c_ptr), intent(in), value :: app
end subroutine app_clear

subroutine app_create(app, appname, intra_comm, root) bind(c,name="app_create")
    import :: c_ptr,c_int,c_char
    type(c_ptr), intent(in), value :: app
    character(kind=c_char), dimension(*), intent(in) :: appname
    type(c_ptr), intent(in), value :: intra_comm
    integer(c_int), intent(in), value :: root
end subroutine app_create

function app_add_boundary(app) bind(c,name="app_add_boundary")
    import :: c_ptr
    type(c_ptr), intent(in), value :: app
    type(c_ptr) :: app_add_boundary
end function app_add_boundary

subroutine app_register_field(app,fieldname,ncomp,location,iotype,units) bind(c,name="app_register_field")
    import :: c_ptr,c_char,c_int
    type(c_ptr), intent(in), value :: app
    character(kind=c_char), dimension(*), intent(in) :: fieldname
    integer(c_int), intent(in), value :: ncomp
    integer(kind(NodeCentered)), intent(in), value :: location
    integer(kind(IncomingDofs)), intent(in), value :: iotype
    character(kind=c_char), dimension(*), intent(in) :: units
end subroutine app_register_field

subroutine app_start_coupling(app, inter_comm) bind(c,name="app_start_coupling")
    import :: c_ptr
    type(c_ptr), intent(in), value :: app, inter_comm
end subroutine app_start_coupling

subroutine app_exchange_solu(app, getter, setter, time, user_data) bind(c,name="app_exchange_solu")
    import :: c_ptr,c_funptr,c_double
    type(c_ptr),    intent(in), value :: app
    type(c_funptr), intent(in), value :: getter,setter
    real(c_double), intent(in), value :: time
    type(c_ptr),    intent(in), value :: user_data
end subroutine app_exchange_solu

subroutine app_stop_coupling(app) bind(c,name="app_stop_coupling")
    import :: c_ptr
    type(c_ptr), intent(in), value :: app
end subroutine app_stop_coupling

!-----------------------------------------
! interfaces for communicator
!-----------------------------------------

function cm_socket_new(as_mastrer,np,master_ip,master_port) bind(c,name="cm_socket_new")
    import :: c_ptr,c_int,c_char
    integer(c_int), intent(in), value :: as_mastrer, np
    character(kind=c_char), dimension(*), intent(in) :: master_ip
    integer(c_int), intent(in), value :: master_port
    type(c_ptr) :: cm_socket_new
end function cm_socket_new

function cm_mpi_new(mpi_comm)  bind(c,name="cm_mpi_new")
    import :: c_ptr,c_int
    integer(c_int), intent(in), value :: mpi_comm
    type(c_ptr) :: cm_mpi_new
end function cm_mpi_new

subroutine cm_set_constant(cm,cname,cvalue) bind(c,name="cm_set_constant")
    import :: c_ptr,c_int,c_char
    type(c_ptr), intent(in), value :: cm
    character(kind=c_char), dimension(*), intent(in) :: cname
    integer(c_int), intent(in), value :: cvalue
end subroutine cm_set_constant

subroutine cm_set_pointer(cm,cname,cvalue) bind(c,name="cm_set_pointer")
    import :: c_ptr,c_char
    type(c_ptr), intent(in), value :: cm
    character(kind=c_char), dimension(*), intent(in) :: cname
    type(c_ptr), intent(in), value :: cvalue
end subroutine cm_set_pointer

subroutine cm_set_function(cm,cname,cvalue) bind(c,name="cm_set_function")
    import :: c_ptr,c_char,c_funptr
    type(c_ptr), intent(in), value :: cm
    character(kind=c_char), dimension(*), intent(in) :: cname
    type(c_funptr), intent(in), value :: cvalue
end subroutine cm_set_function

subroutine cm_delete(cm) bind(c,name="cm_delete")
    import :: c_ptr
    type(c_ptr), intent(inout) :: cm
end subroutine cm_delete

function cm_rank(cm) bind(c,name="cm_rank")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: cm
    integer(c_int) :: cm_rank
end function cm_rank

function cm_size(cm) bind(c,name="cm_size")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: cm
    integer(c_int) :: cm_size
end function cm_size

!-----------------------------------------
! interfaces for Boundary
!-----------------------------------------

function bd_new() bind(c,name="bd_new")
    import :: c_ptr
    type(c_ptr) :: bc_new
end function bd_new

subroutine bd_delete(bd) bind(c,name="bd_delete")
    import :: c_ptr
    type(c_ptr), intent(inout) :: bd
end subroutine bd_delete

subroutine bd_clear(bd) bind(c,name="bd_clear")
    import :: c_ptr
    type(c_ptr), intent(in), value :: bd
end subroutine bd_clear

subroutine bd_set_user_id(bd, id) bind(c,name="bd_set_user_id")
    import :: c_ptr,c_int
    type(c_ptr),    intent(in), value :: bd
    integer(c_int), intent(in), value :: id
end subroutine bd_set_user_id

function bd_get_user_id(bd) bind(c,name="bd_get_user_id")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: bd
    integer(c_int)                 :: bd_get_user_id
end function bd_get_user_id

subroutine bd_reserve(bd, max_node, max_face, max_face_nodes) bind(c,name="bd_reserve")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: bd
    integer(c_int), intent(in), value :: max_node,max_face,max_face_nodes
end subroutine bd_reserve

function bd_add_node(bd,x,y,z,unique_id) bind(c,name="bd_add_node")
    import :: c_ptr,c_int,c_long_long,c_double
    type(c_ptr), intent(in), value :: bd
    real(kind=c_double), intent(in), value :: x,y,z
    integer(c_long_long), intent(in), value :: unique_id
    integer(c_int) :: bd_add_node
end function bd_add_node

function bd_add_face(bd,ftype,nnodes,fnodes) bind(c,name="bd_add_node")
    import :: c_ptr,c_int,c_long_long
    type(c_ptr), intent(in), value :: bd
    integer(kind(FT_BAR2)), intent(in), value :: ftype
    integer(c_int), intent(in), value :: nnodes
    integer(c_int), dimension(*), intent(in) :: fnodes
    integer(c_int) :: bd_add_face
end function bd_add_face

subroutine bd_set_face_centroid(bd,face,cx,cy,cz) bind(c,name="bd_set_face_centroid")
    import :: c_ptr,c_int,c_double
    type(c_ptr), intent(in), value :: bd
    integer(c_int), intent(in), value :: face
    real(kind=c_double), value :: cx,cy,cz
end subroutine bd_set_face_centroid

subroutine bd_set_face_area(bd,face,sx,sy,sz) bind(c,name="bd_set_face_area")
    import :: c_ptr,c_int,c_double
    type(c_ptr), intent(in), value :: bd
    integer(c_int), intent(in), value :: face
    real(kind=c_double), value :: sx,sy,sz
end subroutine bd_set_face_area

subroutine bd_set_node_coords(bd,node,x,y,z) bind(c,name="bd_set_node_coords")
    import :: c_ptr,c_int,c_double
    type(c_ptr), intent(in), value :: bd
    integer(c_int), intent(in), value :: node
    real(kind=c_double), intent(in), value :: x,y,z
end subroutine bd_set_node_coords

function bd_face_num(bd) bind(c,name="bd_face_num")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: bd
    integer(c_int) :: bd_face_num
end function bd_face_num

function bd_node_num(bd) bind(c,name="bd_node_num")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: bd
    integer(c_int) :: bd_face_num
end function bd_node_num

subroutine bd_compute_metrics(bd,basied_angle_deg) bind(c,name="bd_compute_metrics")
    import :: c_ptr,c_double
    type(c_ptr), intent(in), value :: bd
    real(kind=c_double), intent(in), value :: basied_angle_deg
end subroutine bd_compute_metrics

function bd_face_normal(bd,face) bind(c,name="bd_face_normal")
    import :: c_ptr,c_int !,c_double
    type(c_ptr), intent(in), value :: bd
    integer(c_int), intent(in), value :: face
    type(c_ptr) :: bd_face_normal
end function bd_face_normal

function bd_face_area(bd,face) bind(c,name="bd_face_area")
    import :: c_ptr,c_int,c_double
    type(c_ptr), intent(in), value :: bd
    integer(c_int), intent(in), value :: face
    real(kind=c_double) :: bd_face_area
end function bd_face_area

function bd_face_centroid(bd,face) bind(c,name="bd_face_centroid")
    import :: c_ptr,c_int !,c_double
    type(c_ptr), intent(in), value :: bd
    integer(c_int), intent(in), value :: face
    type(c_ptr) :: bd_face_centroid
end function bd_face_centroid

function bd_node_coords(bd,node) bind(c,name="bd_node_coords")
    import :: c_ptr,c_int!,c_double
    type(c_ptr), intent(in), value :: bd
    integer(c_int), intent(in), value :: node
    type(c_ptr) :: bd_node_coords
end function bd_node_coords

! get IndexSet object of boundary nodes.
function bd_nodes(bd) bind(c,name="bd_nodes")
    import :: c_ptr
    type(c_ptr), intent(in), value :: bd
    type(c_ptr) :: bd_nodes
end function bd_nodes

! get MeshConnectivity object of boundary faces.
function bd_face_nodes(bd) bind(c,name="bd_face_nodes")
    import :: c_ptr
    type(c_ptr), intent(in), value :: bd
    type(c_ptr) :: bd_face_nodes
end function bd_face_nodes

! get KdTree object of boundary nodes.
function bd_kdtree(bd) bind(c,name="bd_kdtree")
    import :: c_ptr
    type(c_ptr), intent(in), value :: bd
    type(c_ptr) :: bd_kdtree
end function bd_kdtree

subroutine bd_read_gmsh(bd,filename) bind(c,name="bd_read_gmsh")
    import :: c_ptr,c_char
    type(c_ptr), intent(in), value :: bd
    character(kind=c_char), dimension(*), intent(in) :: filename
end subroutine bd_read_gmsh

!-----------------------------------------
! interfaces for MeshConnectivity
!-----------------------------------------

function mc_new() bind(c,name="mc_new")
    import :: c_ptr
    type(c_ptr) :: mc_new
end function mc_new

subroutine mc_delete(mc) bind(c,name="mc_delete")
    import :: c_ptr
    type(c_ptr), intent(inout) :: mc
end subroutine mc_delete

subroutine mc_clear(mc) bind(c,name="mc_clear")
    import :: c_ptr
    type(c_ptr), intent(in), value :: mc
end subroutine mc_clear

subroutine mc_reserve(mc,max_nrow,max_ndata) bind(c,name="mc_reserve")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: mc
    integer(c_int), intent(in), value :: max_nrow,max_ndata
end subroutine mc_reserve

subroutine mc_push_back(mc,n,pdata) bind(c,name="mc_push_back")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: mc
    integer(c_int), intent(in), value :: n
    integer(c_int), dimension(*), intent(in) :: pdata
end subroutine mc_push_back

function mc_ia(mc) bind(c,name="mc_ia")
    import :: c_ptr!,c_int
    type(c_ptr), intent(in), value :: mc
    type(c_ptr) :: mc_ia
end function mc_ia

function mc_ja(mc) bind(c,name="mc_ja")
    import :: c_ptr!,c_int
    type(c_ptr), intent(in), value :: mc
    type(c_ptr) :: mc_ja
end function mc_ja

function mc_nrow(mc) bind(c,name="mc_nrow")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: mc
    integer(c_int) :: mc_nrow
end function mc_nrow

function mc_ndata(mc) bind(c,name="mc_ndata")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: mc
    integer(c_int) :: mc_ndata
end function mc_ndata

function mc_row_size(mc,row) bind(c,name="mc_row_size")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: mc
    integer(c_int), value :: row
    integer(c_int) :: mc_row_size
end function mc_row_size

function mc_row_data(mc,row,idata) bind(c,name="mc_row_data")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: mc
    integer(c_int), intent(in), value :: row,idata
    integer(c_int) :: mc_row_data
end function mc_row_data

!-----------------------------------------
! interfaces for IndexSet
!-----------------------------------------

function is_new() bind(c,name="is_new")
    import :: c_ptr
    type(c_ptr) :: is_new
end function is_new

subroutine is_delete(is) bind(c,name="is_delete")
    import :: c_ptr
    type(c_ptr), intent(inout) :: is
end subroutine is_delete

subroutine is_clear(is) bind(c,name="is_clear")
    import :: c_ptr
    type(c_ptr), intent(in), value :: is
end subroutine is_clear

function is_contains(is, unique_id) bind(c,name="is_contains")
    import :: c_ptr,c_long_long,c_int
    type(c_ptr), intent(in), value :: is
    integer(c_long_long), intent(in), value :: unique_id
    integer(c_int) :: is_contains
end function is_contains

function is_add(is, unique_id) bind(c,name="is_add")
    import :: c_ptr,c_int,c_long_long
    type(c_ptr), intent(in), value :: is
    integer(c_long_long), intent(in), value :: unique_id
    integer(c_int) :: is_add
end function is_add

function is_l2g(is, ilocal) bind(c,name="is_l2g")
    import :: c_ptr,c_int,c_long_long
    type(c_ptr), intent(in), value :: is
    integer(c_int), intent(in), value :: ilocal
    integer(c_long_long) :: is_l2g
end function is_l2g

function is_g2l(is, iglobal) bind(c,name="is_g2l")
    import :: c_ptr,c_int,c_long_long
    type(c_ptr), intent(in), value :: is
    integer(c_long_long), intent(in), value :: iglobal
    integer(c_int) :: is_g2l
end function is_g2l

function is_size(is) bind(c,name="is_size")
    import :: c_ptr,c_int
    type(c_ptr), intent(in), value :: is
    integer(c_int) :: is_size
end function is_size

function is_glist(is) bind(c,name="is_glist")
    import :: c_ptr,c_int,c_long_long
    type(c_ptr), intent(in), value :: is
    type(c_ptr) :: is_glist
end function is_glist

!-----------------------------------------
! interfaces for KdTree
!-----------------------------------------


!-----------------------------------------
! interfaces for DynamicMatrix
!-----------------------------------------


!-----------------------------------------
! interfaces for DynamicVector
!-----------------------------------------


end interface

contains

subroutine cm_set_function_mpi_initialized(cm,func)
    use iso_c_binding
    type(c_ptr), intent(in), value :: cm
    procedure(mpi_initialized_t), intent(in), pointer :: func
    call cm_set_function(cm,'MPI_Initialized'//char(0),c_funloc(func))
end subroutine cm_set_function_mpi_initialized
subroutine cm_set_function_mpi_init(cm,func)
    use iso_c_binding
    type(c_ptr), intent(in), value :: cm
    procedure(mpi_init_t), intent(in), pointer :: func
    call cm_set_function(cm,'MPI_Init'//char(0),c_funloc(func))
end subroutine cm_set_function_mpi_init
subroutine cm_set_function_mpi_comm_disconnect(cm,func)
    use iso_c_binding
    type(c_ptr), intent(in), value :: cm
    procedure(mpi_comm_disconnect_t), intent(in), pointer :: func
    call cm_set_function(cm,'MPI_Comm_disconnect'//char(0),c_funloc(func))
end subroutine cm_set_function_mpi_comm_disconnect
subroutine cm_set_function_mpi_comm_rank(cm,func)
    use iso_c_binding
    type(c_ptr), intent(in), value :: cm
    procedure(mpi_comm_rank_t), intent(in), pointer :: func
    call cm_set_function(cm,'MPI_Comm_rank'//char(0),c_funloc(func))
end subroutine cm_set_function_mpi_comm_rank
subroutine cm_set_function_mpi_comm_size(cm,func)
    use iso_c_binding
    type(c_ptr), intent(in), value :: cm
    procedure(mpi_comm_size_t), intent(in), pointer :: func
    call cm_set_function(cm,'MPI_Comm_size'//char(0),c_funloc(func))
end subroutine cm_set_function_mpi_comm_size
subroutine cm_set_function_mpi_send(cm,func)
    use iso_c_binding
    type(c_ptr), intent(in), value :: cm
    procedure(mpi_send_t), intent(in), pointer :: func
    call cm_set_function(cm,'MPI_Send'//char(0),c_funloc(func))
end subroutine cm_set_function_mpi_send
subroutine cm_set_function_mpi_recv(cm,func)
    use iso_c_binding
    type(c_ptr), intent(in), value :: cm
    procedure(mpi_recv_t), intent(in), pointer :: func
    call cm_set_function(cm,'MPI_Recv'//char(0),c_funloc(func))
end subroutine cm_set_function_mpi_recv
subroutine app_exchange_solu_wapper(app, getter, setter, time, user_data)
    use iso_c_binding
    type(c_ptr),    intent(in), value :: app
    procedure(get_boundary_field_t), intent(in), pointer :: getter
    procedure(set_boundary_field_t), intent(in), pointer :: setter
    real(c_double), intent(in), value :: time
    type(c_ptr),    intent(in), value :: user_data
    call app_exchange_solu(app,c_funloc(getter),c_funloc(setter),time,user_data)
end subroutine app_exchange_solu_wapper
!#endif

end module easyfsi