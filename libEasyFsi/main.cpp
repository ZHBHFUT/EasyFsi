#include <iostream>
#include <numbers>

//#include "SocketCommunicator.hpp"
#include "Interpolator.hpp"
#include "Boundary.hpp"
#include "Field.hpp"

using namespace EasyLib;

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv)
{
    //std::cout << "Hello World!" << std::endl;

    //int np = 1;
    //const char* name = "unnamed";
    //bool master = false;
    //
    //if (argc <= 1)return -1;

    //for (int i = 1; i < argc; ++i) {
    //    if      (strcmp(argv[i], "-np") == 0) {
    //        ++i;
    //        np = atoi(argv[i]);
    //    }
    //    else if (strcmp(argv[i], "-name") == 0) {
    //        ++i;
    //        name = argv[i];
    //    }
    //    else if (strcmp(argv[i], "-master") == 0) {
    //        master = true;
    //    }
    //}
    
    //EasyLib::SocketCommunicator comm;
    ////comm.init(name, master, np);
    //comm.init(argc, (const char**)argv);
    //
    //std::cout << "np = " << comm.size() << ", rank = " << comm.rank() << '\n';
    //std::cout.flush();
    //
    //std::string msg_send = "Hello";
    //std::string msg_recv;
    //
    //for (int i = 0; i < comm.size(); ++i) {
    //    if (i == comm.rank())continue;
    //    comm.send(msg_send, i, 0);
    //    std::cout << "send \""<< msg_send << "\" to [" << i << "]\n";
    //}
    //
    //for (int i = 0; i < comm.size(); ++i) {
    //    if (i == comm.rank())continue;
    //    comm.recv(msg_recv, i, 0);
    //    std::cout << "recv from [" << i << "] " << msg_recv << '\n';
    //    std::cout.flush();
    //}
    //
    //comm.disconnect();

    Boundary src, des;

    src.read_gmsh("quad_s.msh");
    des.read_gmsh("quad_f.msh");
    src.set_name("solid");
    des.set_name("fluid");

    Boundary::write_to_file("solid.bounds", 1, &src);
    Boundary::write_to_file("fluid.bounds", 1, &des);

    Interpolator interp;
    interp.add_source_boundary(src);
    interp.add_target_boundary(des);

    //interp.compute_interp_coeff();
    interp.compute_interp_coeff(LocalXPS);
    //interp.read_coefficients("test.coeff.dat");

    FieldInfo finfo_s, finfo_f;
    finfo_s.name = "Field";
    finfo_s.units = "";
    finfo_s.ncomp = 1;
    finfo_s.location = NodeCentered;
    finfo_s.iotype = OutgoingDofs;

    finfo_f.name = "Field";
    finfo_f.units = "";
    finfo_f.ncomp = 1;
    finfo_f.location = NodeCentered;
    finfo_f.iotype = IncomingDofs;

    Field field_s, field_f;
    field_s.info = &finfo_s;
    field_f.info = &finfo_f;

    field_s.data.resize(src.node_num(), field_s.info->ncomp);
    field_f.data.resize(des.node_num(), field_f.info->ncomp);

    auto dist = []([[maybe_unused]]double x, [[maybe_unused]] double y, double z) {
        constexpr double pi = std::numbers::pi;
        //double r = sqrt(x * x + y * y + z * z);
        //double phi = atan2(y, x);
        //double theta = atan2(z, r);
        //
        return 0.1 * sin(pi * x) * sin(pi * x) * sin(pi * y) * sin(pi * y) + z;
        //return 0.1 * sin(phi) * sin(theta);
        //return z * 0.1;
    };

    const auto nns = src.node_num();
    //std::vector<double> fields_s(2 * nns);
    
    auto& xyzs = src.node_coords();
    for (int i = 0; i < nns; ++i) {
        field_s.data(i, 0) = dist(xyzs[i].x, xyzs[i].y, xyzs[i].z);
    }

    Field* const ps[1] = { &field_s };
    Field* pf[1] = { &field_f };
    std::span<Field* const> fs(ps, 1);
    std::span< Field*> ff(pf, 1);
    interp.interp_dofs_s2t(fs, ff);

    auto& xyzf = des.node_coords();
    std::vector<double> uf_p(des.node_num());
    double err = 0, delt;
    for (int i = 0; i < des.node_num(); ++i) {
        uf_p[i] = dist(xyzf[i].x, xyzf[i].y, xyzf[i].z);
        delt = field_f.data(i, 0) - uf_p[i];
        err += delt * delt;
        //fields_f[2 * i + 1] = log10(fabs(delt) + 1E-20);
    }
    err = sqrt(err) / (des.node_num());
    std::cout << "ERR = " << err << '\n';

    const char* names[] = {"Field"};
    const double* a[] = { field_s.data.data() };
    const double* b[] = { field_f.data.data() };
    interp.write_tecplot("test.plt", 1, names, a, b);
    interp.save_coefficients("test.coeff.dat");

    return 0;
}
