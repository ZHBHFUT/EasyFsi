#include "EasyFsi.h"

static Application* app;
static Communicator*    inter_comm;
static MPICommunicator* intra_comm;
static Boundary* bd0;

void init()
{
    inter_comm = cm_socket_new(); // 创建socket通信子
    intra_comm = cm_mpi_new(MPI_COMM_WORLD); // 创建应用程序内部多个进程间的通信子
    app = app_new("cfd",intra_comm,0);
    
    // 定义边界
    bd0 = app_add_boundary(app);
    // bd_add_node(bd0,x1,y1,z1,id1);
    // bd_add_node(bd0,x2,y2,z2,id2);
    // ...
    // bd_add_node(bd0,xN,yN,zN,idN);
    // bd_add_face(bd,ft1,nn1,nodes1);
    // bd_add_face(bd,ft2,nn2,nodes2);
    // ...
    // bd_add_face(bd,ftM,nnM,nodesM);
    
    // 定义耦合物理量
    app_register_field(app,"displacement",3,NodeCentered,OutgoingDofs,"m");
    app_register_field(app,"force",3,NodeCentered,IncomingLoads,"N");
    
    // 开始耦合
    app_start_coupling(app,inter_comm);
}

void get_boundary_field(const Boundary* bd, const char* name, int ncomp, FieldLocation loc, double* data, void* user_data)
{
    // 如何获取边界数据
    // ...
}

void set_boundary_field(const Boundary* bd, const char* name, int ncomp, FieldLocation loc, const double* data, void* user_data)
{
    // 如何设置边界值
    // ...
}

// 应用程序的求解函数
void solve()
{
    // ...
    
    // 交换数据
    app_exchange_solu(app,get_boundary_field,set_boundary_field,time,nullptr);
    
    //...
}
