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
//! @file       pybind.cpp
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#ifdef _DEBUG
#define Py_DEBUG
#endif
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include "../libEasyFsi/Application.hpp"
#include "../libEasyFsi/MPICommunicator.hpp"
#include "../libEasyFsi/SocketCommunicator.hpp"

//PY_VERSION

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
        .def("__bool__", [](const type& a) { return !a.empty(); }, "Check whether the vector is non-empty")
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
        .def("to_buffer", to_buffer)
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
        .def("inverse", mat_inv)
        .def("apply", mat_dot_vec)
        .def("apply_add", mat_dot_add_vec)
        .def("apply", mat_dot_mat)
        .def("apply_add", mat_dot_add_mat)
        .def("__bool__", [](const type& a) { return !a.empty(); }, "Check whether the array is non-empty")
        .def("__getitem__", get_item)
        .def("__setitem__", set_item)
        .def("__repr__", print, "Return the canonical string representation of this array.");
    ;
}

void bind_fieldinfo(py::module m)
{
    using type = EasyLib::FieldInfo;

    py::class_<type>(m, "FieldInfo")
        .def(py::init([](py::str name, py::str units, int ncomp, EasyLib::FieldLocation location, EasyLib::FieldIO iotype) { return new type{ (std::string)name, (std::string)units, ncomp, location, iotype }; }))
        .def_readonly("name", &type::name)
        .def_readonly("units", &type::units)
        .def_readonly("ncomp", &type::ncomp)
        .def_readonly("location", &type::location)
        .def_readonly("iotype", &type::iotype)
        ;
}
void bind_field(py::module m)
{
    using type = EasyLib::Field;
    using value_type = typename decltype(type::data)::value_type;
    auto to_buffer = [](type& v) {
        return py::buffer_info{
            v.data.data(),                               // Pointer to buffer
            sizeof(value_type),                          // Size of one scalar
            py::format_descriptor<value_type>::format(), // Python struct-style format descriptor
            2,                                           // Number of dimensions
            {v.data.nrow(), v.data.ncol()},              // Buffer dimensions
            {sizeof(value_type) * v.data.ncol(), sizeof(value_type)} // Strides (in bytes) for each index
        };
    };
    auto get_item = [](const type& a, py::tuple idx) {
        if (idx.size() != 2)throw py::index_error();
        auto i = idx.begin()->cast<int>();
        auto j = (idx.begin() + 1)->cast<int>();
        if (i < 0 || i >= a.data.extent(0) || j < 0 || j >= a.data.extent(1))
            throw py::index_error();
        return a.data(i, j);
    };
    auto set_item = [](type& a, py::tuple idx, double value) {
        if (idx.size() != 2)throw py::index_error();
        auto i = idx.begin()->cast<int>();
        auto j = (idx.begin() + 1)->cast<int>();
        if (i < 0 || i >= a.data.extent(0) || j < 0 || j >= a.data.extent(1))
            throw py::index_error();
        a.set(i, j, value);
    };
    auto fill = [](type& a, double value) { a.data.fill(value); };
    auto get_name = [](type& a) { return a.info ? a.info->name : std::string{}; };
    auto get_units = [](type& a) { return a.info ? a.info->units : std::string{}; };
    auto get_ncomp = [](type& a) { return a.info ? a.info->ncomp : 0; };
    auto get_iotype = [](type& a) { return a.info ? a.info->iotype : EasyLib::OutgoingDofs; };
    auto get_loc = [](type& a) { return a.info ? a.info->location : EasyLib::NodeCentered; };
    auto from_buffer = [](type& a, py::buffer b) {
        if (a.info && !a.info->is_outgoing())
            throw std::runtime_error("Field is not outgoing!");

        // Request a buffer descriptor from Python
        auto info = b.request();

        // Some basic validation checks ...
        if (info.format != py::format_descriptor<value_type>::format())
            throw std::runtime_error("Incompatible format: expected a double array!");
        if (info.ndim != 2)
            throw std::runtime_error("Incompatible buffer dimension!");
        if (info.strides[1] != sizeof(double))
            throw std::runtime_error("Incompatible stride!");
        if (info.shape[0] != a.data.nrow() || info.shape[1] != a.data.ncol())
            throw std::runtime_error("Incompatible shape!");

        std::memcpy(a.data.data(), info.ptr, sizeof(value_type) * a.data.numel());
    };

    py::class_<type>(m, "Field", py::buffer_protocol())
        .def_buffer(to_buffer)
        .def(py::init<>())
        .def("fill", fill)
        .def("set_values", from_buffer, "Set outgoing field values from py::buffer")
        .def("get_values", to_buffer, "Get field values as py::buffer_info")
        .def("__getitem__", get_item)
        .def("__setitem__", set_item)
        .def_property_readonly("name", get_name)
        .def_property_readonly("units", get_units)
        .def_property_readonly("ncomp", get_ncomp)
        .def_property_readonly("iotype", get_iotype)
        .def_property_readonly("location", get_loc)
        .def("__bool__", [](const type& a) { return !a.data.empty(); }, "Check whether the field is non-empty")
        ;
}

class PyCommunicator : public EasyLib::Communicator
{
public:
    using Communicator::Communicator; // Inherit constructors

    //void init(int argc, const char** argv) override { PYBIND11_OVERRIDE_PURE(void, EasyLib::Communicator, init, argc, argv); }

    void set_constant(const char* name, int   value  )override { PYBIND11_OVERRIDE_PURE(void, EasyLib::Communicator, set_constant, name, value); }
    void set_constant(const char* name, void* pointer)override { PYBIND11_OVERRIDE_PURE(void, EasyLib::Communicator, set_constant, name, pointer); }
    void set_function(const char* name, void* func   )override { PYBIND11_OVERRIDE_PURE(void, EasyLib::Communicator, set_constant, name, func); }

    int rank()const noexcept override { PYBIND11_OVERRIDE_PURE(int, EasyLib::Communicator, rank); }

    int size()const noexcept override { PYBIND11_OVERRIDE_PURE(int, EasyLib::Communicator, size); }

    void send(const void* data, int count, EasyLib::DataType type, int dest_rank, int tag)override { PYBIND11_OVERRIDE_PURE(void, EasyLib::Communicator, send, data, count, type, dest_rank, tag); }

    void recv(void* data, int count, EasyLib::DataType type, int src_rank, int tag)override { PYBIND11_OVERRIDE_PURE(void, EasyLib::Communicator, recv, data, count, type, src_rank, tag); }

    void disconnect() override { PYBIND11_OVERRIDE_PURE(void, EasyLib::Communicator, disconnect); }
};

struct MatView
{
    int nrow{ 0 }, ncol{ 0 };
    double* data{ nullptr };
    bool readonly{ false };
};

void bind_matview(py::module m)
{
    using type = MatView;
    auto to_buffer = [](type& mat) {
        return py::buffer_info{
            mat.data,                                     // Pointer to buffer
            sizeof(double),                               // Size of one scalar
            py::format_descriptor<double>::format(),      // Python struct-style format descriptor
            2,                                            // Number of dimensions
            { mat.nrow, mat.ncol },                       // Buffer dimensions
            { sizeof(double) * mat.ncol, sizeof(double) }, // Strides (in bytes) for each index
            mat.readonly
        };
    };
    auto get_item = [](const type& a, py::tuple idx) {
        if (idx.size() == 1) {
            auto i = idx.begin()->cast<int>();
            if (i < 0 || i >= a.nrow * a.ncol)throw py::index_error();
            return a.data[i];
        }
        else if (idx.size() == 2) {
            auto i = idx.begin()->cast<int>();
            auto j = (idx.begin() + 1)->cast<int>();
            if (i < 0 || i >= a.nrow || j < 0 || j >= a.ncol)
                throw py::index_error();
            return a.data[i * a.ncol + j];
        }
        throw py::index_error();
    };
    auto set_item = [](type& a, py::tuple idx, double value) {
        if (a.readonly)throw std::runtime_error("matrix is read-only!");
        if (idx.size() == 1) {
            auto i = idx.begin()->cast<int>();
            if(i<0 || i>=a.nrow*a.ncol)throw py::index_error();
            a.data[i] = value;
        }
        else if (idx.size() == 2) {
            auto i = idx.begin()->cast<int>();
            auto j = (idx.begin() + 1)->cast<int>();
            if (i < 0 || i >= a.nrow || j < 0 || j >= a.ncol)
                throw py::index_error();
            a.data[i * a.ncol + j] = value;
        }
        else
            throw py::index_error();
    };
    py::class_<type>(m, "MatView", py::buffer_protocol())
        .def_buffer(to_buffer)
        .def(py::init<>())
        .def("fill", [](type& m, double value) { std::fill(m.data, m.data + m.nrow * m.ncol, value); })
        .def("__getitem__", get_item)
        .def("__setitem__", set_item)
        .def_readonly("nrow", &type::nrow)
        .def_readonly("ncol", &type::ncol)
        ;
}

static std::unordered_map<const EasyLib::Application*, py::function> py_get_funcs;
static std::unordered_map<const EasyLib::Application*, py::function> py_set_funcs;
void py_get_bound_field(const EasyLib::Application* app, const EasyLib::Boundary* bd, const char* name, int ncomp, EasyLib::FieldLocation loc, double* data, void* user_data)
{
    auto it = py_get_funcs.find(app);
    if (it != py_get_funcs.end()) {
        int nrow = loc == EasyLib::NodeCentered ? bd->nnode() : bd->nface();
        // def py_get_bound_field_func(app, bound, field_name, location, mat, user_data)
        it->second(app, bd, name, loc, MatView{ nrow, ncomp, data, false }, (PyObject*)user_data);
    }
    else {
        throw std::runtime_error("field getter function is missing!");
    }
}
void py_set_bound_field(const EasyLib::Application* app, const EasyLib::Boundary* bd, const char* name, int ncomp, EasyLib::FieldLocation loc, const double* data, void* user_data)
{
    auto it = py_get_funcs.find(app);
    if (it != py_get_funcs.end()) {
        // def py_get_bound_field_func(app, bound, field_name, location, mat, user_data)
        int nrow = loc == EasyLib::NodeCentered ? bd->nnode() : bd->nface();
        it->second(app, bd, name, ncomp, loc, MatView{ nrow, ncomp, const_cast<double*>(data), true }, (PyObject*)user_data);
    }
    else {
        throw std::runtime_error("field setter function is missing!");
    }
}

PYBIND11_MODULE(EasyFsi, m) {
    using namespace EasyLib;

    m.doc() = "pybind11 example plug-in"; // optional module doc string

    //m.def("add", &add, "A function that adds two numbers");

    //m.def("set_func",  [](py::function f) {g_func = f; });
    //m.def("test_func", [](int i, int j) {return g_func(i, j); });

    //--- export constants

    py::enum_<EasyLib::FieldLocation>(m, "FieldLocation")
        .value("NodeCentered", EasyLib::FieldLocation::NodeCentered, "field stored at node")
        .value("CellCentered", EasyLib::FieldLocation::CellCentered, "field stored at cell")
        .export_values();
    py::enum_<EasyLib::FieldIO>(m, "FieldIO")
        .value("IncomingDofs",  EasyLib::FieldIO::IncomingDofs,  "field is incoming DOF")
        .value("IncomingLoads", EasyLib::FieldIO::IncomingLoads, "field is incoming load")
        .value("OutgoingDofs",  EasyLib::FieldIO::OutgoingDofs,  "field is outgoing DOF")
        .value("OutgoingLoads", EasyLib::FieldIO::OutgoingLoads, "field is outgoing load")
        .export_values();
    py::enum_<EasyLib::FaceTopo>(m, "FaceTopo")
        .value("BAR2"   , EasyLib::FaceTopo::BAR2   , "two-node line element")
        .value("BAR3"   , EasyLib::FaceTopo::BAR3   , "three-node line element")
        .value("TRI3"   , EasyLib::FaceTopo::TRI3   , "three-node triangle element")
        .value("TRI6"   , EasyLib::FaceTopo::TRI6   , "six-node triangle element")
        .value("QUAD4"  , EasyLib::FaceTopo::QUAD4  , "four-node quadrilateral element")
        .value("QUAD8"  , EasyLib::FaceTopo::QUAD8  , "eight-node quadrilateral element")
        .value("POLYGON", EasyLib::FaceTopo::POLYGON, "general polygon element")
        .export_values();
    py::enum_<EasyLib::ZoneTopo>(m, "ZoneTopo")
        .value("POINTS"  , EasyLib::ZoneTopo::ZT_POINTS , "points cloud")
        .value("CURVE"   , EasyLib::ZoneTopo::ZT_CURVE  , "curve")
        .value("SURFACE" , EasyLib::ZoneTopo::ZT_SURFACE, "surface")
        .value("VOLUME"  , EasyLib::ZoneTopo::ZT_VOLUME , "3D volume")
        .export_values();
    py::enum_<EasyLib::ZoneShape>(m, "ZoneShape")
        .value("POINT"    , EasyLib::ZoneShape::ZS_POINT   , "zone is a single point")
        .value("COLINEAR" , EasyLib::ZoneShape::ZS_COLINEAR, "zone is a straight line")
        .value("COPLANER" , EasyLib::ZoneShape::ZS_COPLANER, "zone is a plane")
        .value("GENERAL"  , EasyLib::ZoneShape::ZS_GENERAL , "zone is a general 3d surface")
        .export_values();
    //InterpolationMethod
    py::enum_<EasyLib::InterpolationMethod>(m, "InterpolationMethod")
        .value("LocalXPS"  , EasyLib::InterpolationMethod::LocalXPS,   "using local spline method")
        .value("Projection", EasyLib::InterpolationMethod::Projection, "using projection method, elements must exist for source boundary")
        .value("Mapping"   , EasyLib::InterpolationMethod::Mapping,    "using geometric mapping method, elements must exist for all boundary")
        .value("Automatic" , EasyLib::InterpolationMethod::Automatic,  "select interpolation method automatically")
        .export_values();

    // dynamic array
    bind_dynamic_array1<int_l>(m, "VectorInt");
    bind_dynamic_array1<int_g>(m, "VectorIntG");

    // vector and matrix
    bind_vector(m);
    bind_matrix(m);

    // MatView
    bind_matview(m);

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

    // Field
    bind_field(m);

    // FieldInfo
    bind_fieldinfo(m);

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
    auto bd_add_face = [](Boundary& bd, FaceTopo type, py::tuple nodes) {
        std::vector<int_l> list; list.reserve(nodes.size());
        for (auto x : nodes)list.push_back(x.cast<int_l>());
        return bd.add_face(type, (int)list.size(), list.data());
    };
    auto bd_get_field = [](Boundary& bd, py::str name) { return &bd.field(name); };
    py::class_<Boundary>(m, "Boundary")
        .def(py::init<>())
        .def("clear", &Boundary::clear)
        .def("reserve", &Boundary::reserve)
        //.def("set_name", &Boundary::set_name)
        .def_property("name", &Boundary::name, &Boundary::set_name)
        .def_property("user_id", &Boundary::user_id, &Boundary::set_user_id)
        .def("add_node", bd_add_node)
        .def("add_face", bd_add_face)
        .def("set_face_cent", &Boundary::set_face_cent)
        .def("set_face_area", &Boundary::set_face_area)
        .def("compute_metics", &Boundary::compute_metics)
        .def_property_readonly("kdtree", &Boundary::kdtree, py::return_value_policy::reference)
        //.def("compute_global_xps_matrix", &Boundary::compute_global_xps_matrix)
        //.def("compute_global_xps_interp_coeff", &Boundary::compute_global_xps_interp_coeff)
        //.def("compute_local_xps_interp_coeff", &Boundary::compute_local_xps_interp_coeff)
        //.def("compute_project_interp_coeff", &Boundary::compute_project_interp_coeff)
        .def("load", &Boundary::load)
        .def("save", &Boundary::save)
        .def_property_readonly("topo", &Boundary::topo)
        .def_property_readonly("shape", &Boundary::shape)
        .def("contains_polygon", &Boundary::contains_polygon)
        .def("contains_high_order_face", &Boundary::contains_high_order_face)
        .def("all_high_order", &Boundary::all_high_order)
        .def_property_readonly("nnode", &Boundary::nnode)
        .def_property_readonly("nface", &Boundary::nface)
        .def_property_readonly("nelem", &Boundary::nelem)
        .def("get_field", bd_get_field, py::return_value_policy::reference)
        .def("node_coords",   [](Boundary& bound, int_l node) { return &bound.node_coords().at(node); }, py::return_value_policy::reference)
        .def("face_centroid", [](Boundary& bound, int_l face) { return &bound.face_centroids().at(face); }, py::return_value_policy::reference)
        .def("face_area",     [](Boundary& bound, int_l face) { return bound.face_areas().at(face); })
        .def("face_normal",   [](Boundary& bound, int_l face) { return &bound.face_normals().at(face); }, py::return_value_policy::reference)
        .def("face_type",     [](Boundary& bound, int_l face) { return bound.face_types().at(face); })
        .def("register_field", &Boundary::register_field)
        .def("remove_all_field", &Boundary::remove_all_field)
        ;

    auto comm_send_buffer = [](Communicator& comm, py::buffer b, int dest_rank, int tag) {
        // Request a buffer descriptor from Python
        auto info = b.request();

        if (info.size > std::numeric_limits<int>::max())
            throw std::overflow_error("Communicator::send(), length overflow!");

        int count = static_cast<int>(info.size);

        // send 
        if      (info.format == py::format_descriptor<int8_t>::format())
            comm.send(reinterpret_cast<const int8_t*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<int16_t>::format())
            comm.send(reinterpret_cast<const int16_t*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<int32_t>::format())
            comm.send(reinterpret_cast<const int32_t*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<int64_t>::format())
            comm.send(reinterpret_cast<const int64_t*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<uint8_t>::format())
            comm.send(reinterpret_cast<const uint8_t*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<uint16_t>::format())
            comm.send(reinterpret_cast<const uint16_t*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<uint32_t>::format())
            comm.send(reinterpret_cast<const uint32_t*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<uint64_t>::format())
            comm.send(reinterpret_cast<const uint64_t*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<double>::format())
            comm.send(reinterpret_cast<const double*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<float>::format())
            comm.send(reinterpret_cast<const float*>(info.ptr), count, dest_rank, tag);
        else if (info.format == py::format_descriptor<char>::format())
            comm.send(reinterpret_cast<const char*>(info.ptr), count, dest_rank, tag);
        else
            throw std::runtime_error("Communicator::send(), unsupported data type!");
    };
    auto comm_recv_buffer = [](Communicator& comm, py::buffer b, int src_rank, int tag) {
        // Request a buffer descriptor from Python
        auto info = b.request();
        if (info.readonly)throw std::runtime_error("Communicator::recv(), dest storage is read-only!");

        if (info.size > std::numeric_limits<int>::max())
            throw std::overflow_error("Communicator::send(), length overflow!");

        int count = static_cast<int>(info.size);

        // recv 
        if      (info.format == py::format_descriptor<int8_t>::format())
            comm.recv(reinterpret_cast<int8_t*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<int16_t>::format())
            comm.recv(reinterpret_cast<int16_t*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<int32_t>::format())
            comm.recv(reinterpret_cast<int32_t*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<int64_t>::format())
            comm.recv(reinterpret_cast<int64_t*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<uint8_t>::format())
            comm.recv(reinterpret_cast<uint8_t*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<uint16_t>::format())
            comm.recv(reinterpret_cast<uint16_t*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<uint32_t>::format())
            comm.recv(reinterpret_cast<uint32_t*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<uint64_t>::format())
            comm.recv(reinterpret_cast<uint64_t*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<double>::format())
            comm.recv(reinterpret_cast<double*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<float>::format())
            comm.recv(reinterpret_cast<float*>(info.ptr), count, src_rank, tag);
        else if (info.format == py::format_descriptor<char>::format())
            comm.recv(reinterpret_cast<char*>(info.ptr), count, src_rank, tag);
        else
            throw std::runtime_error("Communicator::recv(), unsupported data type!");
    };
    py::class_<Communicator, PyCommunicator>(m, "Communicator")
        .def(py::init<>())
        .def("set_constant", [](Communicator& comm, const char* name, int   value) {comm.set_constant(name, value); })
        .def("set_constant", [](Communicator& comm, const char* name, void* value) {comm.set_constant(name, value); })
        .def("set_function", [](Communicator& comm, const char* name, void* value) {comm.set_function(name, value); })
        .def("rank", &Communicator::rank)
        .def("size", &Communicator::size)
        .def("disconnect", &Communicator::disconnect)
        .def("send", comm_send_buffer)
        .def("recv", comm_recv_buffer)
        .def("send", [](Communicator& comm, const IndexSet& v, int dest, int tag) {comm.send(v, dest, tag); })
        .def("send", [](Communicator& comm, const DynamicVector& v, int dest, int tag) {comm.send(v, dest, tag); })
        .def("send", [](Communicator& comm, const DynamicMatrix& v, int dest, int tag) {comm.send(v, dest, tag); })
        .def("send", [](Communicator& comm, const MeshConnectivity& v, int dest, int tag) {comm.send(v, dest, tag); })
        .def("send", [](Communicator& comm, const Boundary& v, int dest, int tag) {comm.send(v, dest, tag); })
        .def("send", [](Communicator& comm, const std::string& v, int dest, int tag) {comm.send(v, dest, tag); })
        .def("recv", [](Communicator& comm, IndexSet& v, int src, int tag) {comm.recv(v, src, tag); })
        .def("recv", [](Communicator& comm, DynamicVector& v, int src, int tag) {comm.recv(v, src, tag); })
        .def("recv", [](Communicator& comm, DynamicMatrix& v, int src, int tag) {comm.recv(v, src, tag); })
        .def("recv", [](Communicator& comm, MeshConnectivity& v, int src, int tag) {comm.recv(v, src, tag); })
        .def("recv", [](Communicator& comm, Boundary& v, int src, int tag) {comm.recv(v, src, tag); })
        .def("recv", [](Communicator& comm, std::string& v, int src, int tag) {comm.recv(v, src, tag); })
        ;

    // MPI Communicator
    auto mm_set_const = [](MPICommunicator& comm, const char* name, int value) { comm.set_constant(name, value); };
    py::class_<MPICommunicator, Communicator>(m, "MPICommunicator")
        .def(py::init<>())
        .def(py::init<int, int, int>())
        .def_property("MPI_DATATYPE_NULL",   [](const MPICommunicator& c) {return c.get_constant("MPI_DATATYPE_NULL"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_DATATYPE_NULL", value); })
        .def_property("MPI_INT16",     [](const MPICommunicator& c) {return c.get_constant("MPI_INT16_T"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_INT16_T", value); })
        .def_property("MPI_INT32",     [](const MPICommunicator& c) {return c.get_constant("MPI_INT32_T"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_INT32_T", value); })
        .def_property("MPI_INT64",     [](const MPICommunicator& c) {return c.get_constant("MPI_INT64_T"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_INT64_T", value); })
        .def_property("MPI_FLOAT",     [](const MPICommunicator& c) {return c.get_constant("MPI_FLOAT"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_FLOAT", value); })
        .def_property("MPI_DOUBLE",    [](const MPICommunicator& c) {return c.get_constant("MPI_DOUBLE"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_DOUBLE", value); })
        .def_property("MPI_CHAR",      [](const MPICommunicator& c) {return c.get_constant("MPI_CHAR"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_CHAR", value); })
        //.def_property("MPI_SUCCESS",   [](const MPICommunicator& c) {return c.get_constant("MPI_SUCCESS"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_SUCCESS", value); })
        //.def_property("MPI_COMM_WORLD",[](const MPICommunicator& c) {return c.get_constant("MPI_COMM_WORLD"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_COMM_WORLD", value); })
        //.def_property("MPI_COMM_SELF", [](const MPICommunicator& c) {return c.get_constant("MPI_COMM_SELF"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_COMM_SELF", value); })
        //.def_property("MPI_COMM_NULL", [](const MPICommunicator& c) {return c.get_constant("MPI_COMM_NULL"); }, [](MPICommunicator& c, int value) { c.set_constant("MPI_COMM_NULL", value); })
        ;

    // SocketComm
    py::class_<SocketCommunicator, Communicator>(m, "SocketCommunicator")
        .def(py::init<>())
        .def("init", [](SocketCommunicator& sc, bool as_root, int np, py::str master_ip, unsigned short port, int timeout_sec) { std::string s = master_ip; sc.init(as_root, np, s.c_str(), port, timeout_sec); })
        ;

    // Application
    auto set_app_field_funcs = [](Application& app, py::function get_field, py::function set_field) {
        py_get_funcs.emplace(&app, get_field);
        py_set_funcs.emplace(&app, set_field);
        app.set_field_function(&py_get_bound_field, &py_set_bound_field);
    };
    py::class_<Application>(m, "Application")
        .def(py::init<>())
        .def(py::init<const char*>())
        .def(py::init<const char*, Communicator&, int>())
        .def("clear", &Application::clear)
        .def("add_coupled_boundary", [](Application& app) {return &app.add_coupled_boundary(); }, py::return_value_policy::reference)
        .def("boundary_num", &Application::boundary_num)
        .def("boundary", [](Application& app, int b) {return app.boundary(b); }, py::return_value_policy::reference)
        .def("set_field_function", set_app_field_funcs)
        .def("register_field", &Application::register_field)
        .def("start_coupling", &Application::start_coupling)
        .def("exchange_solution", [](Application& app, double time, py::object obj) {app.exchange_solution(time, obj.ptr()); })
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
        .def("compute_interp_coeff", [](Interpolator& it, InterpolationMethod method) { it.compute_interp_coeff(method); })
        .def("compute_interp_coeff", [](Interpolator& it) { it.compute_interp_coeff(); })
        .def("save_coefficients", &Interpolator::save_coefficients)
        .def("load_coefficients", &Interpolator::load_coefficients)
        .def("interp_all_dofs_s2t", &Interpolator::interp_all_dofs_s2t)
        .def("interp_all_load_t2s", &Interpolator::interp_all_load_t2s)
        .def("interp_dofs_s2t", [](Interpolator& it, py::str name) {it.interp_dofs_s2t(std::string(name).c_str()); })
        .def("interp_load_t2s", [](Interpolator& it, py::str name) {it.interp_load_t2s(std::string(name).c_str()); })
        .def("interp_modal_results", [](Interpolator& it, py::str ifile, py::str ofile) {it.interp_modal_results(std::string(ifile).c_str(), std::string(ofile).c_str()); })
        ;
}
