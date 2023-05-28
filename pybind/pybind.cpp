#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include "../libEasyFsi/Application.hpp"
#include "../libEasyFsi/MPICommunicator.hpp"
#include "../libEasyFsi/SocketCommunicator.hpp"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(EasyLib::DynamicArray<EasyLib::int_l, 1>);
PYBIND11_MAKE_OPAQUE(EasyLib::DynamicArray<EasyLib::int_g, 1>);
PYBIND11_MAKE_OPAQUE(EasyLib::DynamicVector);
PYBIND11_MAKE_OPAQUE(EasyLib::DynamicMatrix);

template<typename T>
void bind_dynamic_array1(py::module m, const char* name)
{
    //using namespace EasyLib;
    using type = EasyLib::DynamicArray<T, 1>;
    auto get_item = [](const type& a, int i) {
        if (i < 0 && i >= a.size())throw py::index_error();
        return a[i];
    };
    auto set_item = [](type& a, int i, T value) {
        if (i < 0 && i >= a.size())throw py::index_error();
        a[i] = value;
    };
    auto print = [](const type& a) {
        std::ostringstream s;
        s << '{';
        for (int i = 0; i < a.size(); ++i)
            if (i == 0)
                s << a[i];
            else
                s << ',' << a[i];
        s << '}';
        return s.str();
    };
    py::class_<type>(m, name)
        .def(py::init<>())
        .def(py::init([](int n) {return type(n); }))
        .def(py::init([](int n, T value) {return type(n, value); }))
        .def("rank", &type::rank)
        .def("numel", &type::numel)
        .def("empty", &type::empty)
        .def("size", &type::size)
        .def("clear", &type::clear)
        .def("reserve", &type::reserve)
        .def("resize", [](type& a, size_t new_size) {a.resize(new_size); })
        .def("fill", &type::fill)
        .def("swap", &type::swap)
        .def("swap_elements", &type::swap_elements)
        .def("copy_elements", &type::copy_elements)
        .def("push_back", [](type& a, T value) {a.push_back(value); })
        .def("pop_back", &type::pop_back)
        //.def("front", &type::front)
        //.def("back", &type::back)
        .def("__bool__", [](const type& a) { return !a.empty(); }, "Check whether the array is nonempty")
        .def("__len__", &type::size)
        .def("__getitem__", get_item)
        .def("__setitem__", set_item)
        .def("__repr__", print, "Return the canonical string representation of this array.");
        ;
}
template<typename T>
void bind_dynamic_array2(py::module m, const char* name)
{
    using type = EasyLib::DynamicArray<T, 2>;
    auto get_item = [](const type& a, py::tuple idx) {
        if (idx.size() != 2)throw py::index_error();
        auto i = idx.begin()->cast<int>();
        auto j = (idx.begin() + 1)->cast<int>();
        if (i < 0 || i >= a.extent(0) || j < 0 || j >= a.extent(1))
            throw py::index_error();
        return a(i, j);
    };
    auto set_item = [](type& a, py::tuple idx, T value) {
        if (idx.size() != 2)throw py::index_error();
        auto i = idx.begin()->cast<int>();
        auto j = (idx.begin() + 1)->cast<int>();
        if (i < 0 || i >= a.extent(0) || j < 0 || j >= a.extent(1))
            throw py::index_error();
        a(i, j) = value;
    };
    auto print = [](const type& a) {
        std::ostringstream s;
        s << '{';
        for (int i = 0; i < a.extent(0); ++i) {
            if (i > 0)s << ",\n";
            s << '{';
            for (int j = 0; j < a.extent(1); ++j) {
                if (j == 0)
                    s << a(i, j);
                else
                    s << ',' << a(i, j);
            }
            s << "}";
        }
        s << '}';
        return s.str();
    };
    py::class_<type>(m, name)
        .def(py::init<>())
        .def(py::init([](int n) {return type(n); }))
        .def(py::init([](int n, T value) {return type(n, value); }))
        .def("rank", &type::rank)
        .def("numel", &type::numel)
        .def("empty", &type::empty)
        .def("clear", &type::clear)
        .def("reserve", &type::reserve)
        .def("resize", &type::resize)
        .def("fill", &type::fill)
        .def("swap", &type::swap)
        .def("swap_elements", &type::swap_elements)
        .def("copy_elements", &type::copy_elements)
        .def("__bool__", [](const type& a) { return !a.empty(); }, "Check whether the array is nonempty")
        .def("__getitem__", get_item)
        .def("__setitem__", set_item)
        .def("__repr__", print, "Return the canonical string representation of this array.");
    ;
}

template<typename T>
py::buffer_info span_to_buffer(std::span<T> s)
{
    return py::buffer_info{
        s.data(),                           // Pointer to buffer
        sizeof(T),                          // Size of one scalar
        py::format_descriptor<T>::format(), // Python struct-style format descriptor
        1,                                  // Number of dimensions
        {s.size()},                         // Buffer dimensions
        {sizeof(T)}                         // Strides (in bytes) for each index
    };
}
template<typename T>
std::span<T> buffer_to_span(py::buffer b)
{
    // Request a buffer descriptor from Python
    auto info = b.request();

    // Some basic validation checks ...
    if (info.format != py::format_descriptor<T>::format())
        throw std::runtime_error("Incompatible element type!");
    if (info.ndim != 1)
        throw std::runtime_error("Incompatible buffer dimension!");
    if (info.strides[0] != sizeof(T))
        throw std::runtime_error("Incompatible stride!");

    // make span
    return std::span<T>(reinterpret_cast<T*>(info.ptr), (size_t)info.shape[0]);
}
template<typename T>
auto get_span_item(std::span<T> s, size_t i)
{
    if (i < 0 && i >= s.size())throw py::index_error();
    return s[i];
}
template<typename T>
void set_span_item(std::span<T> s, size_t i, T value)
{
    if (i < 0 && i >= s.size())throw py::index_error();
    s[i] = value;
}
template<typename T>
std::string print_span(std::span<T> a)
{
    std::ostringstream s;
    s << '{';
    for (int i = 0; i < a.size(); ++i)
        if (i == 0)
            s << a[i];
        else
            s << ',' << a[i];
    s << '}';
    return s.str();
}
void bind_span(py::module m, const char* name)
{
    // std::span<int_l>
    {
        using value_type = EasyLib::int_l;
        using type = std::span<value_type>;
        py::class_<type>(m, "SpanInt", py::buffer_protocol())
            .def_buffer(&span_to_buffer<value_type>)
            .def(py::init<>())
            .def(py::init(&buffer_to_span<value_type>))
            .def("size", &type::size)
            .def("__bool__", [](const type& a) { return !a.empty(); }, "Check whether the span is nonempty")
            .def("__len__", &type::size)
            .def("__getitem__", get_span_item<value_type>)
            .def("__setitem__", set_span_item<value_type>)
            .def("__repr__", print_span<value_type>, "Return the canonical string representation of this vector.");
            ;
    }

}

void bind_vector(py::module m)
{
    //using namespace EasyLib;
    using type = EasyLib::DynamicVector;
    auto get_item = [](const type& a, int i) {
        if (i < 0 && i >= a.size())throw py::index_error();
        return a[i];
    };
    auto set_item = [](type& a, int i, double value) {
        if (i < 0 && i >= a.size())throw py::index_error();
        a[i] = value;
    };
    auto print = [](const type& a) {
        std::ostringstream s;
        s << '{';
        for (int i = 0; i < a.size(); ++i)
            if (i == 0)
                s << a[i];
            else
                s << ',' << a[i];
        s << '}';
        return s.str();
    };
    auto vec_min = [](const type& v) {return v.min(); };
    auto vec_max = [](const type& v) {return v.max(); };
    auto to_buffer = [](type& v) {
        return py::buffer_info{
            v.data(),                                          // Pointer to buffer
            sizeof(type::value_type),                          // Size of one scalar
            py::format_descriptor<type::value_type>::format(), // Python struct-style format descriptor
            1,                                                 // Number of dimensions
            {v.size()},                                        // Buffer dimensions
            {sizeof(type::value_type)}                         // Strides (in bytes) for each index
        };
    };
    auto from_buffer = [](py::buffer b) {
        // Request a buffer descriptor from Python
        auto info = b.request();

        // Some basic validation checks ...
        if (info.format != py::format_descriptor<type::value_type>::format())
            throw std::runtime_error("Incompatible format: expected a double array!");
        if (info.ndim != 1)
            throw std::runtime_error("Incompatible buffer dimension!");
        if (info.strides[0] != sizeof(double))
            throw std::runtime_error("Incompatible stride!");

        type v((type::size_type)info.shape[0]);
        std::memcpy(v.data(), info.ptr, sizeof(type::value_type) * v.size());
        return v;
    };
    py::class_<type>(m, "Vector", py::buffer_protocol())
        .def_buffer(to_buffer)
        .def(py::init<>())
        .def(py::init([](int n) {return type(n); }))
        .def(py::init([](int n, double value) {return type(n, value); }))
        .def(py::init(from_buffer))
        .def("rank", &type::rank)
        .def("numel", &type::numel)
        .def("empty", &type::empty)
        .def("size", &type::size)
        .def("clear", &type::clear)
        .def("reserve", &type::reserve)
        .def("resize", [](type& a, size_t new_size) {a.resize(new_size); })
        .def("fill", &type::fill)
        .def("swap", &type::swap)
        .def("swap_elements", &type::swap_elements)
        .def("push_back", [](type& a, double value) {a.push_back(value); })
        .def("pop_back", &type::pop_back)
        //.def("front", &type::front)
        //.def("back", &type::back)
        .def("copy_elements", &type::copy_elements)
        .def("zero", &type::zero)
        .def("norm_sq", &type::norm_sq)
        .def("norm", &type::norm)
        .def("min", vec_min)
        .def("max", vec_max)
        .def("mean", &type::norm)
        .def("dot", &type::dot)
        .def("__bool__", [](const type& a) { return !a.empty(); }, "Check whether the vector is nonempty")
        .def("__len__", &type::size)
        .def("__getitem__", get_item)
        .def("__setitem__", set_item)
        .def("__repr__", print, "Return the canonical string representation of this vector.");
    ;
}
void bind_matrix(py::module m)
{
    using type = EasyLib::DynamicMatrix;
    using vec = EasyLib::DynamicVector;
    auto get_item = [](const type& a, py::tuple idx) {
        if (idx.size() != 2)throw py::index_error();
        auto i = idx.begin()->cast<int>();
        auto j = (idx.begin() + 1)->cast<int>();
        if (i < 0 || i >= a.extent(0) || j < 0 || j >= a.extent(1))
            throw py::index_error();
        return a(i, j);
    };
    auto set_item = [](type& a, py::tuple idx, double value) {
        if (idx.size() != 2)throw py::index_error();
        auto i = idx.begin()->cast<int>();
        auto j = (idx.begin() + 1)->cast<int>();
        if (i < 0 || i >= a.extent(0) || j < 0 || j >= a.extent(1))
            throw py::index_error();
        a(i, j) = value;
    };
    auto print = [](const type& a) {
        std::ostringstream s;
        s << '{';
        for (int i = 0; i < a.extent(0); ++i) {
            if (i > 0)s << ",\n";
            s << '{';
            for (int j = 0; j < a.extent(1); ++j) {
                if (j == 0)
                    s << a(i, j);
                else
                    s << ',' << a(i, j);
            }
            s << "}";
        }
        s << '}';
        return s.str();
    };
    auto mat_inv = [](type& m) { return m.inverse(); };
    auto mat_dot_vec = [](type& m, vec& x, vec& y) { return m.apply(x, y); };
    auto mat_dot_add_vec = [](type& m, vec& x, vec& y) { return m.apply_add(x, y); };
    auto mat_dot_mat = [](type& m, type& x, type& y) { return m.apply(x, y); };
    auto mat_dot_add_mat = [](type& m, type& x, type& y) { return m.apply_add(x, y); };
    auto to_buffer = [](type& mat) {
        return py::buffer_info{
            mat.data(),                                        // Pointer to buffer
            sizeof(type::value_type),                          // Size of one scalar
            py::format_descriptor<type::value_type>::format(), // Python struct-style format descriptor
            2,                                                 // Number of dimensions
            {mat.nrow(), mat.ncol()},                          // Buffer dimensions
            {sizeof(type::value_type) * mat.ncol(), sizeof(type::value_type)} // Strides (in bytes) for each index
        };
    };
    auto from_buffer = [](py::buffer b) {
        // Request a buffer descriptor from Python
        auto info = b.request();

        // Some basic validation checks ...
        if (info.format != py::format_descriptor<type::value_type>::format())
            throw std::runtime_error("Incompatible format: expected a double array!");
        if (info.ndim != 2)
            throw std::runtime_error("Incompatible buffer dimension!");
        if (info.strides[1] != sizeof(double))
            throw std::runtime_error("Incompatible stride!");

        type mat((type::size_type)info.shape[0], (type::size_type)info.shape[1]);
        std::memcpy(mat.data(), info.ptr, sizeof(type::value_type) * mat.numel());
        return mat;
    };
    py::class_<type>(m, "Matrix", py::buffer_protocol())
        .def_buffer(to_buffer)
        .def(py::init<>())
        .def(py::init([](int m, int n) {return type(m, n); }))
        .def(py::init([](int m, int n, double value) {return type(m, n, value); }))
        .def(py::init(from_buffer))
        .def("rank", &type::rank)
        .def("numel", &type::numel)
        .def("empty", &type::empty)
        .def("clear", &type::clear)
        .def("reserve", &type::reserve)
        .def("resize", &type::resize)
        .def("fill", &type::fill)
        .def("swap", &type::swap)
        .def("swap_elements", &type::swap_elements)
        .def("is_square", &type::is_square)
        .def("nrow", &type::nrow)
        .def("ncol", &type::ncol)
        .def("copy_elements", &type::copy_elements)
        .def("zero", &type::zero)
        .def("identity", &type::identity)
        //.def("zero", &type::zero)
        .def("inverse", mat_inv)
        .def("apply", mat_dot_vec)
        .def("apply_add", mat_dot_add_vec)
        .def("apply", mat_dot_mat)
        .def("apply_add", mat_dot_add_mat)
        .def("__bool__", [](const type& a) { return !a.empty(); }, "Check whether the array is nonempty")
        .def("__getitem__", get_item)
        .def("__setitem__", set_item)
        .def("__repr__", print, "Return the canonical string representation of this array.");
    ;
}

//class PyMPIComm : public EasyLib::MPICommunicator
//{
//public:
//
//private:
//
//};

int __stdcall py_MPI_Send(const void* buffer, int count, int datatype, int dest, int tag, int comm)
{
    //
    auto args = py::make_tuple(buffer, count, datatype, dest, tag, comm);
    //TODO:
    // pyfunc(*args)

    return 0;
}
int __stdcall py_MPI_Recv(void* buffer, int count, int datatype, int source, int tag, int comm)
{
    py::make_tuple(buffer, count, datatype, source, tag, comm);
    //TODO:
    return 0;
}

PYBIND11_MODULE(EasyFsi, m) {
    using namespace EasyLib;

    m.doc() = "pybind11 example plug-in"; // optional module doc string

    //m.def("add", &add, "A function that adds two numbers");

    //--- export constants

    py::enum_<EasyLib::FieldLocation>(m, "FieldLocation")
        .value("NodeCentered", EasyLib::FieldLocation::NodeCentered, "field stored at node")
        .value("FaceCentered", EasyLib::FieldLocation::FaceCentered, "field stored at face")
        .value("CellCentered", EasyLib::FieldLocation::CellCentered, "field stored at cell")
        .export_values();
    py::enum_<EasyLib::FieldIO>(m, "FieldIO")
        .value("IncomingDofs ", EasyLib::FieldIO::IncomingDofs,  "field is incoming DOF")
        .value("IncomingLoads", EasyLib::FieldIO::IncomingLoads, "field is incoming load")
        .value("OutgoingDofs ", EasyLib::FieldIO::OutgoingDofs,  "field is outgoing DOF")
        .value("OutgoingLoads", EasyLib::FieldIO::OutgoingLoads, "field is outgoing load")
        .export_values();
    py::enum_<EasyLib::ElementShape>(m, "FaceTopo")
        .value("BAR2 "  , EasyLib::ElementShape::BAR2   , "two-node line element")
        .value("BAR3"   , EasyLib::ElementShape::BAR3   , "three-node line element")
        .value("TRI3 "  , EasyLib::ElementShape::TRI3   , "three-node triangle element")
        .value("TRI6"   , EasyLib::ElementShape::TRI6   , "six-node triangle element")
        .value("QUAD4"  , EasyLib::ElementShape::QUAD4  , "four-node quadrilateral element")
        .value("QUAD8"  , EasyLib::ElementShape::QUAD8  , "eight-node quadrilateral element")
        .value("POLYGON", EasyLib::ElementShape::POLYGON, "general polygon element")
        .export_values();
    py::enum_<EasyLib::ZoneTopo>(m, "ZoneTopo")
        .value("POINTS " , EasyLib::ZoneTopo::ZT_POINTS , "points cloud")
        .value("CURVE"   , EasyLib::ZoneTopo::ZT_CURVE  , "curve")
        .value("SURFACE ", EasyLib::ZoneTopo::ZT_SURFACE, "surface")
        .value("VOLUME"  , EasyLib::ZoneTopo::ZT_VOLUME , "3D volume")
        .export_values();
    py::enum_<EasyLib::ZoneShape>(m, "ZoneShape")
        .value("POINT "   , EasyLib::ZoneShape::ZS_POINT   , "zone is a single point")
        .value("COLINEAR" , EasyLib::ZoneShape::ZS_COLINEAR, "zone is a straight line")
        .value("COPLANER ", EasyLib::ZoneShape::ZS_COPLANER, "zone is a plane")
        .value("GENERAL"  , EasyLib::ZoneShape::ZS_GENERAL , "zone is a general 3d surface")
        .export_values();

    // dynamic array
    bind_dynamic_array1<int_l>(m, "VectorInt");
    bind_dynamic_array1<int_g>(m, "VectorIntG");

    // vector and matrix
    bind_vector(m);
    bind_matrix(m);

    // TinyVector
    using vec3 = TinyVector<double, 3>;
    auto print_vec3 = [](const vec3& v) { return "{x=" + std::to_string(v.x) + ",y=" + std::to_string(v.y) + ",z=" + std::to_string(v.x) + "}"; };
    py::class_<vec3>(m, "Vec3")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def("assign", &vec3::assign)
        .def("fill", &vec3::fill)
        .def("swap", &vec3::swap)
        .def_readwrite("x", &vec3::x)
        .def_readwrite("y", &vec3::y)
        .def_readwrite("z", &vec3::z)
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(py::self *= py::float_())
        .def(py::self /= py::float_())
        .def(py::float_() * py::self)
        .def(py::self * py::float_())
        .def(py::self / py::float_())
        .def(-py::self)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("norm", &vec3::norm)
        .def("norm_sq", &vec3::norm_sq)
        .def("normalize", &vec3::normalize)
        .def("dot", &vec3::dot)
        .def("cross", &vec3::cross)
        .def("distance", &vec3::distance)
        .def("distance_sq", &vec3::distance_sq)
        .def("__repr__", print_vec3)
        ;

    // IndexSet
    auto list2is = [](py::list list) {
        IndexSet::lvec l2g; l2g.reserve(list.size());
        for (auto x : list)l2g.push_back(x.cast<int_g>());
        return IndexSet(std::move(l2g));
    };
    auto is_create = [](IndexSet& is, py::list list) {
        IndexSet::lvec l2g; l2g.reserve(list.size());
        for (auto x : list)l2g.push_back(x.cast<int_g>());
        is.create(std::move(l2g));
    };
    auto is_get_item = [](const IndexSet& is, int_l i) -> int_g {
        if (i < 0 && i >= is.size())throw py::index_error();
        return is[i];
    };
    auto is_to_bool = [](const IndexSet& is) -> bool { return !is.empty(); };
    auto print_is = [](const IndexSet& is) {
        std::ostringstream s;
        s << '{';
        for (int_l i = 0; i < is.size(); ++i)
            if (i == 0)
                s << is[i];
            else
                s << ',' << is[i];
        s << '}';
        return s.str();
    };
    py::class_<EasyLib::IndexSet>(m, "IndexSet")
        .def(py::init<>()) // default constructor
        .def(py::init(list2is))
        .def("clear", &IndexSet::clear)
        .def("create", is_create)
        .def("add", &IndexSet::add)
        .def("contains", &IndexSet::contains)
        .def("find", &IndexSet::find)
        .def("empty", &IndexSet::empty)
        .def("size", &IndexSet::size)
        .def("l2g", &IndexSet::l2g)
        .def("g2l", &IndexSet::g2l)
        .def("__len__", &IndexSet::size)
        .def("__getitem__", is_get_item)
        .def("__bool__", is_to_bool, "Check whether the IndexSet is nonempty")
        .def("__repr__", print_is, "Return the canonical string representation of this IndexSet.");
    ;

    // kdtree
    auto kd_to_bool = [](const KDTree<>& k) {return !k.empty(); };
    auto kd_create = [](KDTree<>& k, DynamicMatrix& coords) {
        if (coords.ncol() != 3)throw py::value_error();
        k.create(coords.data(), coords.nrow(), false);
    };
    py::class_<KDTree<>>(m, "KDTree")
        .def(py::init<>())
        .def("clear", &KDTree<>::clear)
        .def("empty", &KDTree<>::empty)
        .def("size", &KDTree<>::size)
        .def("create", kd_create)
        .def("__bool__", kd_to_bool, "Check whether the kdtree is nonempty")
        ;

    // MeshConnectivity
    auto mc_push_back_tp = [](MeshConnectivity& mc, py::tuple tp) {
        std::vector<int_l> list; list.reserve(tp.size());
        for (auto x : tp)list.push_back(x.cast<int_l>());
        return mc.push_back((int_l)list.size(), list.data());
    };
    auto mc_push_back_list = [](MeshConnectivity& mc, py::list lst) {
        std::vector<int_l> list; list.reserve(lst.size());
        for (auto x : lst)list.push_back(x.cast<int_l>());
        return mc.push_back((int_l)list.size(), list.data());
    };
    auto mc_get_item = [](const MeshConnectivity& mc, py::tuple idx) {
        if (idx.size() != 2)throw py::index_error();
        auto i = idx.begin()->cast<int>();
        auto j = (idx.begin() + 1)->cast<int>();
        if (i < 0 || i >= mc.nrow() || j < 0 || j >= mc.ndata(i))
            throw py::index_error();
        return mc(i, j);
    };
    py::class_<MeshConnectivity>(m, "MeshConnectivity")
        .def(py::init<>())
        .def("clear", &MeshConnectivity::clear)
        .def("reserve", &MeshConnectivity::reserve)
        .def("nrow", &MeshConnectivity::nrow)
        .def("ndata", [](const MeshConnectivity& mc) { return mc.ndata(); })
        .def("ndata", [](const MeshConnectivity& mc, int i) { return mc.ndata(i); })
        .def("empty", &MeshConnectivity::empty)
        .def_property_readonly("ia", &MeshConnectivity::ia)
        .def_property_readonly("ja", &MeshConnectivity::ja)
        .def_static("flip", &MeshConnectivity::flip)
        .def("push_back", mc_push_back_tp)
        .def("push_back", mc_push_back_list)
        .def("__getitem__", mc_get_item)
        ;

    // Boundary
    auto bd_add_node = [](Boundary& bd, double x, double y, double z, int_g unique_id) {return bd.add_node(x, y, z, unique_id); };
    auto bd_add_face = [](Boundary& bd, ElementShape type, py::tuple nodes) {
        std::vector<int_l> list; list.reserve(nodes.size());
        for (auto x : nodes)list.push_back(x.cast<int_l>());
        return bd.add_face(type, (int)list.size(), list.data());
    };
    py::class_<Boundary>(m, "Boundary")
        .def(py::init<>())
        .def("clear", &Boundary::clear)
        .def("reserve", &Boundary::reserve)
        //.def("set_name", &Boundary::set_name)
        .def_property("name", &Boundary::name,&Boundary::set_name)
        .def_property("user_id",&Boundary::user_id, &Boundary::set_user_id)
        .def("add_node", bd_add_node)
        .def("add_face", bd_add_face)
        .def("set_face_cent", &Boundary::set_face_cent)
        .def("set_face_area", &Boundary::set_face_area)
        .def("compute_metics", &Boundary::compute_metics)
        .def_property_readonly("kdtree", &Boundary::kdtree, py::return_value_policy::reference_internal)
        //.def("compute_global_xps_matrix", &Boundary::compute_global_xps_matrix)
        //.def("compute_global_xps_interp_coeff", &Boundary::compute_global_xps_interp_coeff)
        //.def("compute_local_xps_interp_coeff", &Boundary::compute_local_xps_interp_coeff)
        //.def("compute_project_interp_coeff", &Boundary::compute_project_interp_coeff)
        .def("read_gmsh", &Boundary::read_gmsh)
        .def("read_f3d_tec", &Boundary::read_f3d_tec)
        .def_property_readonly("topo", &Boundary::topo)
        .def_property_readonly("shape", &Boundary::shape)
        .def("contains_polygon", &Boundary::contains_polygon)
        .def("contains_high_order_face", &Boundary::contains_high_order_face)
        .def("all_high_order", &Boundary::all_high_order)
        .def_property_readonly("nnode", &Boundary::node_num)
        .def_property_readonly("nface", &Boundary::face_num)
        ;

    // MPI Communicator
    auto mm_set_const = [](MPICommunicator& comm, const char* name, int value) { comm.set_constant(name, value); };
    py::class_<MPICommunicator>(m, "MPICommunicator")
        .def(py::init<>())
        .def(py::init<int, int, int>())
        .def_property("MPI_INT16_T",   [](const MPICommunicator& c) {return c.get_constant("MPI_INT16_T"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_INT16_T", value); })
        .def_property("MPI_INT32_T",   [](const MPICommunicator& c) {return c.get_constant("MPI_INT32_T"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_INT32_T", value); })
        .def_property("MPI_INT64_T",   [](const MPICommunicator& c) {return c.get_constant("MPI_INT64_T"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_INT64_T", value); })
        .def_property("MPI_FLOAT",     [](const MPICommunicator& c) {return c.get_constant("MPI_FLOAT"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_FLOAT", value); })
        .def_property("MPI_DOUBLE",    [](const MPICommunicator& c) {return c.get_constant("MPI_DOUBLE"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_DOUBLE", value); })
        .def_property("MPI_CHAR",      [](const MPICommunicator& c) {return c.get_constant("MPI_CHAR"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_CHAR", value); })
        .def_property("MPI_SUCCESS",   [](const MPICommunicator& c) {return c.get_constant("MPI_SUCCESS"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_SUCCESS", value); })
        .def_property("MPI_COMM_WORLD",[](const MPICommunicator& c) {return c.get_constant("MPI_COMM_WORLD"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_COMM_WORLD", value); })
        .def_property("MPI_COMM_SELF", [](const MPICommunicator& c) {return c.get_constant("MPI_COMM_SELF"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_COMM_SELF", value); })
        .def_property("MPI_COMM_NULL", [](const MPICommunicator& c) {return c.get_constant("MPI_COMM_NULL"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_COMM_NULL", value); })
        ;

    // SocketComm
    py::class_<SocketCommunicator>(m, "SocketCommunicator")
        .def(py::init<>())
        .def("init", [](SocketCommunicator& sc, bool as_root, int np, const char* master_ip , unsigned short port, int timeout_sec) { sc.init(as_root, np, master_ip, port, timeout_sec); })
        .def("disconnect", &SocketCommunicator::disconnect)
        ;

    // Application
    py::class_<Application>(m, "Application")
        .def(py::init<>())
        .def(py::init<const char*>())
        .def(py::init<const char*, Communicator&, int>())
        .def("clear", &Application::clear)
        .def("add_coupled_boundary", &Application::add_coupled_boundary, py::return_value_policy::reference_internal)
        .def("boundary_num", &Application::boundary_num)
        .def("boundary", [](Application& app, int b) {return *app.boundary(b); }, py::return_value_policy::reference_internal)
        .def("register_field", &Application::register_field)
        .def("start_coupling", &Application::start_coupling)
        .def("exchange_solution", &Application::exchange_solution)
        .def("stop_coupling", &Application::stop_coupling)
        .def("save_tecplot", &Application::save_tecplot)
        ;

    // interpolator
    py::class_<Interpolator>(m, "Interpolator")
        .def(py::init<>())
        .def("clear", &Interpolator::clear)
        .def("add_source_boundary", &Interpolator::add_source_boundary)
        .def("add_target_boundary", &Interpolator::add_target_boundary)
        .def("compute_interp_coeff", &Interpolator::compute_interp_coeff)
        .def("save_coefficients", &Interpolator::save_coefficients)
        .def("read_coefficients", &Interpolator::read_coefficients)
        ;
}
