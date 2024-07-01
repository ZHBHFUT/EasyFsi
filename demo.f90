#include "mpif.h"
program demo
    implicit none
    use iso_c_binding,    only : c_int,c_char,c_double,c_long_long,c_ptr,c_funptr,c_funloc,c_null_ptr
    use kinddefs,         only : dp,i8
    
    type(c_ptr), target :: socket_comm
    type(c_ptr), target :: mpi_comm
    type(c_ptr), target :: this_app
    
    logical :: master = .false.
    integer :: n_participator = 2
    character(80) :: ip = "127.0.0.1"
    integer :: port = 40000
    
    integer :: comm_rank,comm_size,ierror
    integer :: nbound, ib
    integer :: it,niter
    double precision :: time
    
    !-------------------------
    ! preprocess 
    !-------------------------
    
    !...
    
    !-------------------------
    ! init FSI
    !-------------------------
    
    ! create socket communicator for inter-communication
    socket_comm = cm_socket_new(int(master,c_int),n_participator,&
                                trim(ip)//char(0),&
                                int(port,c_int))
    
    ! create mpi communicator for intra-communication
    mpi_comm_rank(MPI_COMM_WORLD,comm_rank,ierror)
    mpi_comm_size(MPI_COMM_WORLD,comm_size,ierror)
    mpi_comm = cm_mpi_new(MPI_COMM_WORLD,comm_rank,comm_size)
    ! set constants of MPI
    call cm_set_constant(mpi_comm,'MPI_DATATYPE_NULL'//char(0),   MPI_DATATYPE_NULL)
    call cm_set_constant(mpi_comm,'MPI_INT16_T'//char(0),   MPI_INTEGER2)
    call cm_set_constant(mpi_comm,'MPI_INT32_T'//char(0),   MPI_INTEGER4)
    call cm_set_constant(mpi_comm,'MPI_INT64_T'//char(0),   MPI_INTEGER8)
    call cm_set_constant(mpi_comm,'MPI_FLOAT'//char(0),     MPI_FLOAT)
    call cm_set_constant(mpi_comm,'MPI_DOUBLE'//char(0),    MPI_DOUBLE)
    call cm_set_constant(mpi_comm,'MPI_CHAR'//char(0),      MPI_CHAR)
    
    ! set functions of MPI
    call cm_set_function_mpi_send(mpi_comm,mpi_send_imp)
    call cm_set_function_mpi_recv(mpi_comm,mpi_recv_imp)
    
    ! create application
    this_app = app_new('flow'//char(0), mpi_comm, 0)
    
    ! create fsi boundary
    do ib=1:nbound
        bd = app_add_boundary(this_app)
        
        ! TODO: create bd
        !  bd_add_node(bd,x,y,z,node_g) ! add nodes to boundary
        !...
        !  bd_add_face(bd,TRI3,3,nodes) ! add faces to boundary
        !...
        
        call bd_set_user_id(bd,ib)
        call bd_compute_metrics(bd, 5.0_dp);
    end do
    
    ! register fields
    call app_register_field(this_app,'displacement'//char(0),3,NodeCentered,IncomingDofs, 'Meter'//char(0))
    call app_register_field(this_app,'temperature'//char(0), 1,NodeCentered,IncomingDofs, 'Kelvin'//char(0))
    call app_register_field(this_app,'force'//char(0),       3,NodeCentered,OutgoingLoads,'Newton'//char(0))
    call app_register_field(this_app,'heat'//char(0),        1,NodeCentered,OutgoingLoads,'Watt'//char(0))
    
    ! setup getter and setter function
    call app_set_field_func_wapper(this_app, get_boundary_field, set_boundary_field)
    
    ! start coupling
    call app_start_coupling(this_app, socket_comm)

    !-------------------------
    ! iterate
    !-------------------------
    
    do it = 1:niter
        ! exchange solution
        call app_exchange_solu(this_app, real(time,c_double),c_null_ptr)
        
        ! iterate
        call iterate()
    end do 
    
    !-------------------------
    ! postprocess
    !-------------------------
    
    ! stop coupling and clear FSI
    call app_stop_coupling(this_app)
    if (c_associated(this_app   ))call app_delete(c_loc(this_app))
    if (c_associated(socket_comm))call cm_delete(c_loc(socket_comm))
    if (c_associated(mpi_comm   ))call cm_delete(c_loc(mpi_comm))
    this_app    = c_null_ptr
    socket_comm = c_null_ptr
    mpi_comm    = c_null_ptr
    
    ! ...
    
end program demo

subroutine get_boundary_field(app,bd,fieldname,ncomp,location,data,&
                              user_data) bind(c)
    !DEC$ ATTRIBUTES STDCALL :: get_boundary_field
    use iso_c_binding,           only : c_ptr,c_char,c_int,c_double
    use easyfsi,                 only : bd_get_user_id,NodeCentered
    
    type(c_ptr),                   intent(in), value :: app
    type(c_ptr),                   intent(in), value :: bd
    character(kind=c_char), dimension(*), intent(in) :: fieldname
    integer(kind=c_int),           intent(in), value :: ncomp
    integer(kind(NodeCentered)),   intent(in), value :: location
    real(kind=c_double),   dimension(*), intent(out) :: fielddata
    type(c_ptr),                   intent(in), value :: user_data
    
    integer(c_int) :: ib
    
continue
    ib = bd_get_user_id(bd)  ! the boundary id
    
    ! TODO: compute and fill fielddata from your boundary field
    !do i=1,bc(ib)%nbnode
    !    fielddata(3*(i-1)+1) = bc(ib)%fx(i);
    !    fielddata(3*(i-1)+2) = bc(ib)%fy(i);
    !    fielddata(3*(i-1)+3) = bc(ib)%fz(i);
    !end do
end subroutine

subroutine set_boundary_field(app,bd,fieldname,ncomp,location,fielddata,&
                              user_data) bind(c)
!DEC$ ATTRIBUTES STDCALL :: set_boundary_field
    use iso_c_binding,           only : c_ptr,c_char,c_int,c_double
    use easyfsi,                 only : bd_get_user_id,NodeCentered
    
    type(c_ptr),                   intent(in), value :: app
    type(c_ptr),                   intent(in), value :: bd
    character(c_char),      dimension(*), intent(in) :: fieldname
    integer(c_int),                intent(in), value :: ncomp
    integer(kind(NodeCentered)),   intent(in), value :: location
    real(kind=c_double),    dimension(*), intent(in) :: fielddata
    type(c_ptr),                   intent(in), value :: user_data
    
    integer(c_int) :: ib
    integer        :: i, node_l
    real(dp)       :: x,y,z
    
continue
    ib = bd_get_user_id(bd)  ! the boundary id
    
    ! TODO: update your boundary field with fielddata
    ! 
    !do i=1,bc(ib)%nbnode
    !    bc(ib)%dx(i) = fielddata(3*(i-1)+1);
    !    bc(ib)%dy(i) = fielddata(3*(i-1)+2);
    !    bc(ib)%dz(i) = fielddata(3*(i-1)+3);
    !end do
end subroutine set_boundary_field

function mpi_send_imp(buffer,count,datatype,dest,tag,comm) bind(c)
!DEC$ ATTRIBUTES STDCALL :: mpi_send_imp
    use iso_c_binding, only : c_int,c_ptr
    type(c_ptr), intent(in), value :: buffer
    integer(c_int), intent(in), value :: count
    integer(c_int), intent(in), value :: datatype
    integer(c_int), intent(in), value :: dest
    integer(c_int), intent(in), value :: tag
    integer(c_int), intent(in), value :: comm
    integer(c_int) :: mpi_send_imp
    
    integer ierror
    call mpi_send(buffer,count,datatype,dest,tag,comm,ierror)
    mpi_send_imp = ierror
end function mpi_send_imp
function mpi_recv_imp(buffer,count,datatype,source,tag,comm) bind(c)
!DEC$ ATTRIBUTES STDCALL :: mpi_recv_imp
    use iso_c_binding, only : c_int
    type(c_ptr), intent(in), value :: buffer
    integer(c_int), intent(in), value :: count
    integer(c_int), intent(in), value :: datatype
    integer(c_int), intent(in), value :: source
    integer(c_int), intent(in), value :: tag
    integer(c_int), intent(in), value :: comm
    integer(c_int) :: mpi_recv_imp
    
    integer ierror
    call mpi_recv(buffer,count,datatype,source,tag,comm,MPI_STATUS_IGNORE,ierror)
    mpi_recv_imp = ierror
end function mpi_recv_imp
