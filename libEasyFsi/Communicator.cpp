#include <sstream>

#include "IndexSet.hpp"
#include "MeshConnectivity.hpp"
#include "DynamicVector.hpp"
#include "DynamicMatrix.hpp"
#include "Boundary.hpp"
#include "Field.hpp"
#include "Communicator.hpp"

namespace EasyLib {

    void Communicator::send(const IndexSet& is, int dest_rank, int tag)
    {
        using IS = IndexSet;
        int len = static_cast<int>(is.l2g_.size());
        send(&len, 1, dest_rank, tag); // send length
        if (len > 0)send(is.l2g_.data(), len, dest_rank, tag + 1); // send data
    }

    void Communicator::recv(IndexSet& is, int src_rank, int tag)
    {
        is.clear();

        int len = 0;
        recv(&len, 1, src_rank, tag); // recv length
        if (len > 0) {
            is.l2g_.resize(len);
            recv(is.l2g_.data(), len, src_rank, tag + 1); // recv data
            is.build_g2l_();
        }
    }

    void Communicator::send(const DynamicVector& dv, int dest_rank, int tag)
    {
        int size = static_cast<int>(dv.size());
        send(&size, 1, dest_rank, tag);
        if (size > 0)send(dv.data(), size, dest_rank, tag + 1);
    }

    void Communicator::recv(DynamicVector& dv, int src_rank, int tag)
    {
        dv.clear();
        int size = 0;
        recv(&size, 1, src_rank, tag);
        if (size > 0) {
            dv.resize(size);
            recv(dv.data(), size, src_rank, tag + 1);
        }
    }

    void Communicator::send(const DynamicMatrix& dm, int dest_rank, int tag)
    {
        const int size[2] = { static_cast<int>(dm.extent(0)),static_cast<int>(dm.extent(1)) };
        send(size, 2, dest_rank, tag);
        if (!dm.empty())send(dm.data(), static_cast<int>(dm.numel()), dest_rank, tag + 1);
    }

    void Communicator::recv(DynamicMatrix& dm, int src_rank, int tag)
    {
        using DM = DynamicMatrix;

        int size[2] = { 0 };
        recv(size, 2, src_rank, tag);
        dm.resize(size[0], size[1]);
        if (!dm.empty())recv(dm.data(), static_cast<int_l>(dm.numel()), src_rank, tag + 1);
    }

    void Communicator::send(const MeshConnectivity& mc, int dest_rank, int tag)
    {
        int size[2] = {mc.nrow(), mc.ndata()};
        send(size, 2, dest_rank, tag);
        if (mc.nrow() == 0)return;
        send(mc.ia().data(), mc.nrow() + 1, dest_rank, ++tag);
        if (mc.ndata() == 0)return;
        send(mc.ja().data(), mc.ndata(), dest_rank, ++tag);
    }

    void Communicator::recv(MeshConnectivity& mc, int src_rank, int tag)
    {
        mc.clear();
        int size[2] = { 0 };
        recv(size, 2, src_rank, tag);
        mc.ia_.resize(size[0] + 1, 0);
        if (size[0] <= 1)return;
        recv(mc.ia_.data(), size[1] + 1, src_rank, ++tag);
        mc.ja_.resize(mc.ia_.back());
        if (size[1] > 0)recv(mc.ja_.data(), size[1], src_rank, ++tag);
    }

    void Communicator::send(const Boundary& bound, int dest_rank, int tag)
    {
        ASSERT(bound.mesh_changed_);

        send(bound.name(), dest_rank, tag);
        send(bound.nodes(), dest_rank, ++tag);
        send(bound.face_nodes(), dest_rank, ++tag);

        const auto nn = bound.nodes_.size();
        const auto nf = bound.face_nodes_.nrow();
        const auto nd = Boundary::vec3().size();

        // send node data
        if (bound.node_num() > 0)
            send(bound.node_coords().data()->data(), nd * nn, dest_rank, ++tag);

        // send face data
        if (bound.face_num() > 0) {
            send(bound.face_centroids().data()->data(), nd * nf, dest_rank, ++tag);
            send(bound.face_areas().data(), nf, dest_rank, ++tag);
            send(bound.face_normals().data()->data(), nd * nf, dest_rank, ++tag);
            send(bound.face_types().data(), nf, dest_rank, ++tag);
        }
    }

    void Communicator::recv(Boundary& bound, int src_rank, int tag)
    {
        bound.clear();

        recv(bound.name_, src_rank, tag);
        recv(bound.nodes_, src_rank, ++tag);
        recv(bound.face_nodes_, src_rank, ++tag);
        
        const auto nn = bound.nodes_.size();
        const auto nf = bound.face_nodes_.nrow();
        const auto nd = Boundary::vec3().size();

        // allocate
        bound.node_coords_.resize(3 * nn);
        bound.face_centroids_.resize(3 * nf);
        bound.face_area_.resize(nf);
        bound.face_normal_.resize(3 * nf);
        bound.face_types_.resize(nf);

        // recv node data
        if (nn > 0)
            recv(bound.node_coords_.data()->data(), nd * nn, src_rank, ++tag);

        // recv face data
        if (nf > 0) {
            recv(bound.face_centroids_.data()->data(), nd * nf, src_rank, ++tag);
            recv(bound.face_area_.data(), nf,  src_rank, ++tag);
            recv(bound.face_normal_.data()->data(), nd * nf,  src_rank, ++tag);
            recv(bound.face_types_.data(), nf, src_rank, ++tag);
        }

        bound.mesh_changed_ = true;

        //? Don't compute metrics here.
        //bound.compute_metics();
    }

    void Communicator::send(const std::string& str, int dest_rank, int tag)
    {
        int len = static_cast<int>(str.length());
        send(&len, 1, dest_rank, tag);
        if (!str.empty())send(str.c_str(), len, dest_rank, tag + 1);
    }

    void Communicator::recv(std::string& str, int src_rank, int tag)
    {
        str.clear();
        int len = 0;
        recv(&len, 1, src_rank, tag);
        if (len > 0) {
            str.resize(len);
            recv(str.data(), len, src_rank, tag + 1);
        }
    }

    void Communicator::send(const FieldInfo& info, int dest_rank, int tag)
    {
        std::ostringstream oss;
        oss << info.name << '\n'
            << info.units << '\n'
            << info.ncomp << ' '
            << info.location << ' '
            << info.iotype << ' '
            << info.time << ' '
            << (int)info.is_orphan << ' '
            << (int)info.is_out_of_date
            ;
        return send(oss.str(), dest_rank, tag);
    }

    void Communicator::recv(FieldInfo& info, int src_rank, int tag)
    {
        std::string str;
        recv(str, src_rank, tag);

        std::istringstream iss(str);

        std::getline(iss, info.name);
        std::getline(iss, info.units);
        int loc = NodeCentered, type = IncomingDofs, orphan = 0, out_of_dat = 0;
        iss >> info.ncomp
            >> loc
            >> type
            >> info.time
            >> orphan
            >> out_of_dat;
        info.location       = (FieldLocation)loc;
        info.iotype         = (FieldIO)type;
        info.is_orphan      = orphan;
        info.is_out_of_date = out_of_dat;
    }
}
