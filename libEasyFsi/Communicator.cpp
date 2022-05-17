#include <sstream>

#include "IndexSet.hpp"
#include "MeshConnectivity.hpp"
#include "DynamicVector.hpp"
#include "DynamicMatrix.hpp"
#include "Boundary.hpp"
#include "Field.hpp"
#include "Communicator.hpp"

namespace EasyLib {

    bool Communicator::send(const IndexSet& is, int dest_rank, int tag)
    {
        using IS = IndexSet;
        auto len = static_cast<int_l>(is.l2g_.size());
        if (!this->send(&len, 1, dest_rank, tag))return false; // send length
        return len > 0 ? send(is.l2g_.data(), len, dest_rank, tag + 1) : true; // send data
    }

    bool Communicator::recv(IndexSet& is, int src_rank, int tag)
    {
        using IS = IndexSet;
        is.clear();
        int_l len = 0;
        if (!this->recv(&len, 1, src_rank, tag))return false; // recv length
        if (len <= 0)return true;

        is.l2g_.resize(len);
        if (!recv(is.l2g_.data(), len, src_rank, tag + 1))return false; // recv data

        is.build_g2l_();
        return true;
    }

    bool Communicator::send(const DynamicVector& dv, int dest_rank, int tag)
    {
        using DV = DynamicVector;

        auto size = static_cast<int_l>(dv.size());
        if (!this->send(&size, 1, dest_rank, tag))return false;
        return size > 0 ? send(dv.data(), size, dest_rank, tag + 1) : true;
    }

    bool Communicator::recv(DynamicVector& dv, int src_rank, int tag)
    {
        using DV = DynamicVector;

        dv.clear();
        int_l size = 0;
        if (!this->recv(&size, 1, src_rank, tag))return false;
        dv.resize(size);
        return size > 0 ? recv(dv.data(), size, src_rank, tag + 1) : true;
    }

    bool Communicator::send(const DynamicMatrix& dm, int dest_rank, int tag)
    {
        using DM = DynamicMatrix;

        const int_l size[2] = { static_cast<int_l>(dm.size(0)),static_cast<int_l>(dm.size(1)) };
        if (!this->send(size, 2, dest_rank, tag))return false;
        if (dm.empty())return true;
        return !dm.empty() ? send(dm.data(), static_cast<int_l>(dm.numel()), dest_rank, tag + 1) : true;
    }

    bool Communicator::recv(DynamicMatrix& dm, int src_rank, int tag)
    {
        using DM = DynamicMatrix;

        int_l size[2] = { 0 };
        if (!this->recv(size, 2, src_rank, tag))return false;
        dm.resize(size[0], size[1]);
        return !dm.empty() ? recv(dm.data(), static_cast<int_l>(dm.numel()), src_rank, tag + 1) : true;
    }

    bool Communicator::send(const MeshConnectivity& mc, int dest_rank, int tag)
    {
        using MC = MeshConnectivity;

        int_l size[2] = {mc.nrow(), mc.ndata()};
        if (!this->send(size, 2, dest_rank, tag))return false;
        if (mc.nrow() == 0)return true;
        if (!send(mc.ia().data(), mc.nrow() + 1, dest_rank, ++tag))return false;
        if (mc.ndata() == 0)return true;
        return send(mc.ja().data(), mc.ndata(), dest_rank, ++tag);
    }

    bool Communicator::recv(MeshConnectivity& mc, int src_rank, int tag)
    {
        using MC = MeshConnectivity;
        mc.clear();
        int_l size[2] = { 0 };
        if (!this->recv(size, 2, src_rank, tag))return false;
        mc.ia_.resize(size[0] + 1, 0);
        if (size[0] <= 1)return true;
        if (!recv(mc.ia_.data(), size[1] + 1, src_rank, ++tag))return false;
        mc.ja_.resize(mc.ia_.back());
        return size[1] > 0 ? recv(mc.ja_.data(), size[1], src_rank, ++tag) : true;
    }

    bool Communicator::send(const Boundary& bound, int dest_rank, int tag)
    {
        ASSERT(bound.mesh_changed_);

        send(bound.name(), dest_rank, tag);
        send(bound.nodes(), dest_rank, ++tag);
        send(bound.face_nodes(), dest_rank, ++tag);

        const auto nn = bound.nodes_.size();
        const auto nf = bound.face_nodes_.nrow();
        const auto nd = Boundary::vec3().size();

        // send node data
        if (bound.node_num() > 0) {
            if (!send(bound.node_coords().data()->data(), nd * nn, dest_rank, ++tag))return false;
        }

        // send face data
        if (bound.face_num() > 0) {
            if (!send(bound.face_centroids().data()->data(), nd * nf, dest_rank, ++tag))return false;
            if (!send(bound.face_areas().data(), nf, dest_rank, ++tag))return false;
            if (!send(bound.face_normals().data()->data(), nd * nf, dest_rank, ++tag))return false;
            if (!send(bound.face_types().data(), nf, dest_rank, ++tag))return false;
        }

        return true;
    }

    bool Communicator::recv(Boundary& bound, int src_rank, int tag)
    {
        bound.clear();

        if (!recv(bound.name_, src_rank, tag))return false;
        if (!recv(bound.nodes_, src_rank, ++tag))return false;
        if (!recv(bound.face_nodes_, src_rank, ++tag))return false;
        ++tag;

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
        if (nn > 0) {
            if (!recv(bound.node_coords_.data()->data(), nd * nn, src_rank, ++tag))return false;
        }

        // recv face data
        if (nf > 0) {
            if (!recv(bound.face_centroids_.data()->data(), nd * nf, src_rank, ++tag))return false;
            if (!recv(bound.face_area_.data(), nf,  src_rank, ++tag))return false;
            if (!recv(bound.face_normal_.data()->data(), nd * nf,  src_rank, ++tag))return false;
            if (!recv(bound.face_types_.data(), nf,  src_rank, ++tag))return false;
        }

        bound.mesh_changed_ = true;

        //? Don't compute metrics here.
        //bound.compute_metics();

        return true;
    }

    bool Communicator::send(const std::string& str, int dest_rank, int tag)
    {
        int_l len = static_cast<int_l>(str.length());
        if (!send(&len, 1, dest_rank, tag))return false;
        if (len == 0)return true;
        return send(str.c_str(), len, dest_rank, tag + 1);
    }

    bool Communicator::recv(std::string& str, int src_rank, int tag)
    {
        str.clear();
        int_l len = 0;
        if (!recv(&len, 1, src_rank, tag))return false;
        if (len == 0)return true;
        str.resize(len);
        return recv(str.data(), len, src_rank, tag + 1);
    }

    bool Communicator::send(const FieldInfo& info, int dest_rank, int tag)
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

    bool Communicator::recv(FieldInfo& info, int src_rank, int tag)
    {
        std::string str;
        if (!recv(str, src_rank, tag))return false;

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
        info.location = (FieldLocation)loc;
        info.iotype = (FieldIO)type;
        info.is_orphan = orphan;
        info.is_out_of_date = out_of_dat;

        return true;
    }
}
