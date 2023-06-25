#include "EasyFsi.h"

static Application*     app        = nullptr;
static Communicator*    inter_comm = nullptr;
static MPICommunicator* intra_comm = nullptr;
static Boundary*        bd0        = nullptr;

//-------------------------------------------------------
// define helper functions for field reading and writing
//-------------------------------------------------------

//! @brief function used to reading outgoing fields invoked by Application
void get_boundary_field(const Application* app, const Boundary* bd, const char* name, int ncomp, FieldLocation loc, double* data, void* user_data)
{
    if (loc == NodeCentered){
        int nnode = bd_nnode(bd);
        double* p = data;
        for(int node=0;node<nnode;++node,p+=ncomp){
            // TODO: read nodal field from solver and update output data
            p[0]=???;
            p[1]=???;
            //...
            p[ncomp-1]=???;
        }
    }
    else{
        int nface = bd_nface(bd);
        double* p = data;
        for(int face=0;face<nface;++face,p+=ncomp){
            // TODO: read face field from solver and update output data
            p[0]=???;
            p[1]=???;
            //...
            p[ncomp-1]=???;
        }
    }
}

//! @brief function used to writing incoming fields invoked by Application
void set_boundary_field(const Application* app, const Boundary* bd, const char* name, int ncomp, FieldLocation loc, const double* data, void* user_data)
{
    if (loc == NodeCentered){
        int nnode = bd_nnode(bd);
        double* p = data;
        for(int node=0;node<nnode;++node,p+=ncomp){
            // TODO: update nodal field by input data
            //... = p[0];
            //... = p[1];
            //...
            //... = p[ncomp-1];
        }
    }
    else{
        int nface = bd_nface(bd);
        double* p = data;
        for(int face=0;face<nface;++face,p+=ncomp){
            // TODO: update face field by input data
            //... = p[0];
            //... = p[1];
            //...
            //... = p[ncomp-1];
        }
    }
}

//-------------------------------------------------------
// preprocessing
//-------------------------------------------------------

void init()
{
    // create communicator used to transfer data between different application.
    inter_comm = cm_socket_new(true,2,"127.0.0.1",1234); 
    
    // create MPI communicator used to transfer data between process of this application
    //    arg0   The MPI communicator
    //    arg1   Rank of this process in mpi_comm
    //    arg2   Total process number of mpi_comm
    intra_comm = cm_mpi_new(mpi_comm, rank, size);
    
    // create application
    app = app_new("cfd",intra_comm,0);
    
    // define coupled boundary
    bd0 = app_add_boundary(app);
    // 1) create boundary manually:
    // bd_add_node(bd0,x1,y1,z1,id1);
    // ...
    // bd_add_face(bd,ft1,nn1,nodes1);
    // ...
    //
    // 2) create boundary from Gmsh file:
    // bd_read_gmsh(bd, "???.msh");
    
    bd_compute_metrics(bd0, 5.0);
    
    // define coupled fields:
    app_register_field(app,"displacement",3,NodeCentered,OutgoingDofs ,"m");
    app_register_field(app,"force"       ,3,NodeCentered,IncomingLoads,"N");
    
    // set field reading/writing functions
    app_set_field_func(app,get_boundary_field,set_boundary_field);
    
    // start coupling: get solver infomation
    app_start_coupling(app,inter_comm);
}

//-------------------------------------------------------
// solving
//-------------------------------------------------------

void solve()
{
    // solve this physics one step
    //...
    
    // interpolate and exchange results between applications
    app_exchange_solu(app,time,nullptr);
    
    // other operations
    //...
}

//-------------------------------------------------------
// postprocessing
//-------------------------------------------------------

void post()
{
    app_stop_coupling(app);
    app_delete(&app);
    bd_delete(&bd0);
    cm_delete(&intra_comm);
    cm_delete(&inter_comm);
    
    // other operations
}
