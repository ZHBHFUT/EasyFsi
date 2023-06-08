#pragma once
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
//! @file       TinyVector.hpp
//!             The definition TinyVector class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <cmath> // sqrt, 
#include <initializer_list>
#include <algorithm>
#include <stdexcept>
#include <array>

#include "Assert.hpp"

template<typename T, size_t ND> using StdArray = std::array<T, ND>;

template<typename T, size_t N> class TinyVector : public StdArray<T, N>
{
public:
    using base_type        = StdArray<T, N>;
    using value_type       = typename base_type::value_type;
    using size_type        = typename base_type::size_type;
    using reference        = typename base_type::reference;
    using const_reference  = typename base_type::const_reference;
    using pointer          = T*;
    using const_pointer    = const std::remove_const_t<T>*;

    static_assert(N > 0, "TinyVector<T, N>, vector length must great than zero");
    
    //----------------------------------------
    // constructor
    //----------------------------------------
    
    //! @brief default constructor
    TinyVector() = default;
    //! @brief copy constructor
    TinyVector(const TinyVector&) = default;
    //! assign operator
    TinyVector& operator = (const TinyVector&) = default;

    //! @brief ctor from std::array
    TinyVector(const std::array<value_type, N>& a)noexcept :base_type(a) {}

    //! @brief ctor from c-array
    template<typename T2>
    explicit TinyVector(const T2 (&a)[N])noexcept :base_type()
    {
        auto p = base_type::data();
        for (size_t i = 0; i < N; ++i, ++p)
            (*p) = static_cast<value_type>(a[i]);
    }

    //! @brief ctor from variadic arguments
    template<typename ... Args>
    explicit TinyVector(Args&& ... args)noexcept
        :base_type()
    {
        static_assert(sizeof...(args) == N, "argument number not match the vector length");
        auto it = base_type::begin();
        auto f = [&it](value_type val) {*it = val; ++it; };
        (f(std::forward<Args>(args)), ...);
    }

    //! @brief ctor from initialize list
    template<typename T2>
    TinyVector(std::initializer_list<T2> list)noexcept
        :base_type()
    {
        auto it = list.begin();
        auto p = base_type::data();
        auto q = p + N;
        for (; p != q && it != list.end(); ++it, ++p)
            *p = static_cast<value_type>(*it);
        ASSERT(it == list.end());
        for (; p != q; ++p)*p = value_type{ 0 };
    }

    //----------------------------------------
    // viewer functions
    //----------------------------------------
    
    //! @brief view memory as a TinyVector object.
    static TinyVector& view(pointer data, size_t offset = 0)noexcept
    {
        //static_assert(sizeof(TinyVector) == sizeof(value_type) * N, "TinyVector::view(), type size not agree.");
        return *reinterpret_cast<TinyVector*>(data + offset);
    }

    //! @brief view memory as a TinyVector object.
    static const TinyVector& view(const_pointer data, size_t offset = 0)noexcept
    {
        //static_assert(sizeof(TinyVector) == sizeof(value_type) * N, "TinyVector::view(), type size not agree.");
        return *reinterpret_cast<const TinyVector*>(data + offset);
    }

    //----------------------------------------
    // common methods
    //----------------------------------------

    //! @brief assign vector with values 
    template<typename ... Args>
    void assign(Args&& ... args)noexcept
    {
        static_assert(sizeof...(args) == N, "argument number not match the vector length");
        auto it = base_type::begin();
        auto f = [&it](value_type val) {*it = val; ++it; };
        (f(std::forward<Args>(args)), ...);
    }

    //! @brief fill vector with same value
    void fill(value_type val)noexcept { for (auto& x : base_type)x = val; }

    //! @brief swap to vector
    void swap(TinyVector& b)noexcept { base_type::swap((base_type&)b); }

    reference       max_element()noexcept
    {
        return *std::max_element(base_type::begin(), base_type::end());
    }
    const_reference max_element()const noexcept
    {
        return *std::max_element(base_type::cbegin(), base_type::cend());
    }
    reference       min_element()noexcept
    {
        return *std::min_element(base_type::begin(), base_type::end());
    }
    const_reference min_element()const noexcept
    {
        return *std::min_element(base_type::cbegin(), base_type::cend());
    }

    //----------------------------------------
    // mathematic methods
    //----------------------------------------
    
    //! @brief compute squared L2-norm of the vector
    value_type norm_sq()const noexcept
    {
        value_type l2(0);
        auto p = base_type::data();
        auto q = p + N;
        for (; p != q; ++p)l2 += (*p) * (*p);
        return l2;
    }

    //! @brief compute L2-norm of the vector
    value_type norm()const noexcept
    {
        return std::sqrt(this->norm_sq());
    }

    //! @brief normalize the vector. i.e. make the vector as unit one.
    value_type normalize() noexcept
    {
        value_type l = this->norm_sq();
        if (l != 0) {
            auto p = base_type::data();
            auto q = p + N;
            l = std::sqrt(l);
            for (; p != q; ++p)(*p) /= l;
        }
        return l;
    }

    //! @brief compute dot(a,b) 
    value_type dot(const TinyVector& b)const noexcept
    {
        value_type r{ 0 };
        auto p = base_type::data();
        auto q = p + N;
        auto s = b.data();
        for (; p != q; ++p, ++s)r += (*p) * (*s);
        return r;
    }

    //! Y = a * X + Y 
    TinyVector& add_ax(value_type a, const TinyVector& X)noexcept
    {
        auto p = base_type::data();
        auto q = p + N;
        auto s = X.data();
        for (; p != q; ++p, ++s)(*p) += a * (*s);
        return *this;
    }

    value_type distance_sq(const TinyVector& b)const noexcept
    {
        value_type d{ 0 }, dx;
        auto p = base_type::data();
        auto q = p + N;
        auto s = b.data();
        for (; p != q; ++p, ++s) {
            dx = (*p) - (*s);
            d += dx * dx;
        }
        return d;
    }

    value_type distance(const TinyVector& b)const noexcept
    {
        return std::sqrt(this->distance_sq(b));
    }

    //----------------------------------------
    // assign operator
    //----------------------------------------
    
    //! assign operator by initialize list
    template<typename T2>
    TinyVector& operator = (std::initializer_list<T2> list)
    {
        if (list.size() != N)throw std::out_of_range("TinyVector::operator=(list), length not agree");
        auto it = list.begin();
        auto p  = base_type::data();
        for (size_t i = 0; i < N && it != list.end(); ++i, ++it, ++p)(*p) = static_cast<value_type>(*it);
        return *this;
    }

    //! assign operator by c-array
    template<typename T2>
    TinyVector& operator = (const T2 (&v)[N])noexcept
    {
        auto p = base_type::data();
        for (size_t i = 0; i < N; ++i, ++p)(*p) = static_cast<value_type>(v[i]);
        return *this;
    }

    //! assign operator by std::array
    TinyVector& operator = (const std::array<value_type, N>& v)noexcept
    {
        base_type::operator=(v);
        return *this;
    }

    //! assign operator by value
    //TinyVector& operator = (value_type val)
    //{
    //    auto p = base_type::data();
    //    for (size_t i = 0; i < N; ++i, ++p)(*p) = val;
    //    return *this;
    //}

    //----------------------------------------
    // subscript operator
    //----------------------------------------

    template<size_type I>       value_type& entity()      noexcept { static_assert(I < N, "TinyVector::entity(), index out of range"); return base_type::operator[](I); }
    template<size_type I> const value_type& entity()const noexcept { static_assert(I < N, "TinyVector::entity(), index out of range"); return base_type::operator[](I); }

    //----------------------------------------
    // mathematic operator
    //----------------------------------------

    //! x = scalar * x
    TinyVector& operator *= (const value_type& scale)noexcept
    {
        auto p = base_type::data();
        for (size_t i = 0; i < N; ++i, ++p)(*p) *= scale;
        return *this;
    }

    //! x = x / scalar
    TinyVector& operator /= (const value_type& divisor)noexcept
    {
        auto p = base_type::data();
        for (size_type i = 0; i < N; ++i, ++p)(*p) /= divisor;
        return *this;
    }

    //! x = x + y
    TinyVector& operator += (const TinyVector& y)noexcept
    {
        auto p = base_type::data();
        for (size_type i = 0; i < N; ++i, ++p)(*p) += y[i];
        return *this;
    }

    //! x = x - y
    TinyVector& operator -= (const TinyVector& y)noexcept
    {
        auto p = base_type::data();
        for (size_type i = 0; i < N; ++i, ++p)(*p) -= y[i];
        return *this;
    }

    //! x = x + y
    TinyVector& operator += (const value_type (&y)[N])noexcept
    {
        auto p = base_type::data();
        for (size_type i = 0; i < N; ++i, ++p)(*p) += y[i];
        return *this;
    }
    template<typename VEC>
    TinyVector& operator += (const VEC& y)
    {
        if (y.size() != N) throw std::invalid_argument("TinyVector::operator +=(VEC), vector length not agree");
        auto p = base_type::data();
        for (size_type i = 0; i < N; ++i, ++p) { (*p) += y[i]; }
        return *this;
    }

    //! x = x - y
    TinyVector& operator -= (const value_type(&y)[N])noexcept
    {
        auto p = base_type::data();
        for (size_type i = 0; i < N; ++i, ++p)(*p) -= y[i];
        return *this;
    }
    template<typename VEC>
    TinyVector& operator -= (const VEC& y)
    {
        if (y.size() != N) throw std::invalid_argument("TinyVector::operator -=(VEC), vector length not agree");
        auto p = base_type::data();
        for (size_type i = 0; i < N; ++i, ++p) { (*p) -= y[i]; }
        return *this;
    }

    //! z = x + y
    friend TinyVector operator + (const TinyVector& x, const TinyVector& y)noexcept
    {
        TinyVector z;
        for (size_type i = 0; i < N; ++i)z[i] = x[i] + y[i];
        return z;
    }

    //! z = x - y
    friend TinyVector operator - (const TinyVector& x, const TinyVector& y)noexcept
    {
        TinyVector z;
        for (size_type i = 0; i < N; ++i)z[i] = x[i] - y[i];
        return z;
    }

    //! y = scalar * x
    friend TinyVector operator * (value_type scalar, const TinyVector& x)noexcept
    {
        TinyVector y;
        for (size_type i = 0; i < N; ++i)y[i] = scalar * x[i];
        return y;
    }

    //! y = x * scalar
    friend TinyVector operator * (const TinyVector& x, value_type scalar)noexcept
    {
        TinyVector y;
        for (size_type i = 0; i < N; ++i)y[i] = scalar * x[i];
        return y;
    }

    //! y = x / scalar
    friend TinyVector operator / (const TinyVector& x, value_type scalar)noexcept
    {
        TinyVector y;
        for (size_type i = 0; i < N; ++i)y[i] = x[i] / scalar;
        return y;
    }

    //! y = -x
    friend inline TinyVector operator - (const TinyVector& x)noexcept
    {
        TinyVector y;
        for (size_t i = 0; i < N; ++i)y[i] = -x[i];
        return y;
    }

    //----------------------------------------
    // friend mathematic functions
    //----------------------------------------

    friend value_type dot(const TinyVector& a, const TinyVector& b)noexcept
    {
        auto p = a.data();
        auto q = p + N;
        auto s = b.data();
        value_type ret{ 0 };
        for (; p != q; ++p, ++s)ret += (*p) * (*s);
        return ret;
    }

    friend value_type norm_sq(const TinyVector& a)noexcept
    {
        value_type ret{ 0 };
        auto p = a.data();
        auto q = p + N;
        for (; p != q; ++p)ret += (*p) * (*p);
        return std::sqrt(ret);
    }

    friend value_type norm(const TinyVector& a)noexcept
    {
        return std::sqrt(norm_sq(a));
    }

    friend value_type normalize(TinyVector& a)noexcept
    {
        value_type l = norm_sq(a);
        if (l != 0) {
            auto p = a.data();
            auto q = p + N;
            l = std::sqrt(l);
            for (; p != q; ++p)(*p) /= l;
        }
        return l;
    }

    friend void swap(TinyVector& a, TinyVector& b)noexcept { std::swap((base_type&)a, (base_type&)b); }
};

template<typename T>
struct TinyVector<T, 1>
{
    using value_type             = T;
    using size_type              = size_t;
    using reference              = T&;
    using const_reference        = const std::remove_const_t<T>& ;
    using pointer                = T*;
    using const_pointer          = const std::remove_const_t<T>*;
    using iterator               = typename StdArray<T, 1>::iterator;
    using const_iterator         = typename StdArray<T, 1>::const_iterator;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    value_type x{ 0 };

    //----------------------------------------
    // constructor
    //----------------------------------------
    
    TinyVector() = default;
    TinyVector(const TinyVector&) = default;
    TinyVector& operator = (const TinyVector&) = default;

    explicit TinyVector(value_type X)noexcept :x(X) {}
    explicit TinyVector(const StdArray<T, 1>& X)noexcept :x(X.front()) {}
    explicit TinyVector(value_type (&X)[1])noexcept :x(X[0]) {}

    //----------------------------------------
    // viewer functions
    //----------------------------------------

    //! @brief view memory as a TinyVector object.
    static TinyVector& view(pointer data, size_t offset = 0)noexcept
    {
        return *reinterpret_cast<TinyVector*>(data + offset);
    }

    //! @brief view memory as a TinyVector object.
    static const TinyVector& view(const_pointer data, size_t offset = 0)noexcept
    {
        return *reinterpret_cast<const TinyVector*>(data + offset);
    }

    //----------------------------------------
    // common methods
    //----------------------------------------

    constexpr size_type size()const noexcept { return 1; }

    constexpr size_type max_size() const noexcept { return 1; }

    constexpr bool empty() const noexcept { return false; }

    pointer       data()      noexcept { return &x; }
    const_pointer data()const noexcept { return &x; }

    void assign(value_type X)noexcept { x = X; }

    void fill(value_type val)noexcept { x = val; }

    void swap(TinyVector& b)noexcept { std::swap(x, b.x); }

    reference       max_element()noexcept
    {
        return x;
    }
    const_reference max_element()const noexcept
    {
        return x;
    }
    reference       min_element()noexcept
    {
        return x;
    }
    const_reference min_element()const noexcept
    {
        return x;
    }

    //----------------------------------------
    // mathematic methods
    //----------------------------------------

    value_type norm()const noexcept { return std::abs(x); }

    value_type norm_sq()const noexcept { return x * x; }

    value_type normalize()noexcept
    {
        auto l = std::abs(x);
        x = std::copysign(value_type{ 1 }, x);
        return l;
    }

    value_type dot(const TinyVector& b)const noexcept { return x * b.x; }

    //! Y = a * X + Y
    TinyVector& add_ax(value_type a, const TinyVector& X)noexcept
    {
        x += a * X.x;
        return *this;
    }

    value_type distance_sq(const TinyVector& b)const noexcept
    {
        auto dx = x - b.x;
        return dx * dx;
    }

    value_type distance(const TinyVector& b)const noexcept
    {
        return std::sqrt(this->distance_sq(b));
    }

    //----------------------------------------
    // subscript operator
    //----------------------------------------

    reference       operator[]([[maybe_unused]]size_type i)      noexcept { ASSERT(i == 0); return x; }
    const_reference operator[]([[maybe_unused]]size_type i)const noexcept { ASSERT(i == 0); return x; }

    reference       at(size_type i)
    {
        if (i)throw std::out_of_range("TinyVector<*,1>::at(), index out of range");
        return x;
    }
    const_reference at(size_type i)const
    {
        if (i)throw std::out_of_range("TinyVector<*,1>::at(), index out of range");
        return x;
    }

    template<size_type I>       value_type& entity()      noexcept { static_assert(I == 0, "TinyVector<*,1>::entity(), index out of range"); return x; }
    template<size_type I> const value_type& entity()const noexcept { static_assert(I == 0, "TinyVector<*,1>::entity(), index out of range"); return x; }
    
    //----------------------------------------
    // assign operator
    //----------------------------------------
    
    //! assign operator by value
    //TinyVector& operator = (const value_type& value)noexcept
    //{
    //    x = value;
    //    return *this;
    //}
    
    //! assign operator by value list
    template<typename T2>
    TinyVector& operator = (std::initializer_list<T2> list)
    {
        if (list.size() != size())throw std::out_of_range("TinyVector::operator=(list), length not agree");
        x = static_cast<value_type>(*list.begin());
        return *this;
    }
    //! assign operator by c array
    template<typename T2>
    TinyVector& operator = (const T2(&v)[1])noexcept
    {
        x = static_cast<value_type>(v[0]);
        return *this;
    }
    //! assign operator by std::array
    TinyVector& operator = (const std::array<value_type, 1>& v)noexcept
    {
        x = v[0];
        return *this;
    }
    
    bool operator ==(const TinyVector& b)const noexcept
    {
        return x == b.x;
    }
    bool operator !=(const TinyVector& b)const noexcept
    {
        return x != b.x;
    }

    //----------------------------------------
    // other special operator
    //----------------------------------------

    operator value_type()const noexcept { return x; }

    //----------------------------------------
    // iterator
    //----------------------------------------

    iterator               begin  ()       noexcept { return (iterator(&x, 0)); }
    const_iterator         begin  () const noexcept { return (const_iterator(&x, 0)); }
    iterator               end    ()       noexcept { return (iterator(&x, size())); }
    const_iterator         end    () const noexcept { return (const_iterator(&x, size())); }
    reverse_iterator       rbegin ()       noexcept { return (reverse_iterator(end())); }
    const_reverse_iterator rbegin () const noexcept { return (const_reverse_iterator(end())); }
    reverse_iterator       rend   ()       noexcept { return (reverse_iterator(begin())); }
    const_reverse_iterator rend   () const noexcept { return (const_reverse_iterator(begin())); }
    const_iterator         cbegin () const noexcept { return (begin()); }
    const_iterator         cend   () const noexcept { return (end()); }
    const_reverse_iterator crbegin() const noexcept { return (rbegin()); }
    const_reverse_iterator crend  () const noexcept { return (rend()); }

    //----------------------------------------
    // mathematic operators
    //----------------------------------------
        
    //! x = scalar * x
    TinyVector& operator *= (const value_type& scale)noexcept
    {
        x *= scale;
        return *this;
    }

    //! x = x / scalar
    TinyVector& operator /= (const value_type& divisor)noexcept
    {
        x /= divisor;
        return *this;
    }

    //! x = x + y
    TinyVector& operator += (const TinyVector& Y)noexcept
    {
        x += Y.x;
        return *this;
    }

    //! x = x - y
    TinyVector& operator -= (const TinyVector& Y)noexcept
    {
        x -= Y.x;
        return *this;
    }

    //! z = x + y
    friend TinyVector operator + (const TinyVector& X, const TinyVector& Y)noexcept
    {
        return TinyVector(X.x + Y.x);
    }

    //! z = x - y
    friend TinyVector operator - (const TinyVector& X, const TinyVector& Y)noexcept
    {
        return TinyVector(X.x - Y.x);
    }

    //! y = scalar * x
    friend TinyVector operator * (value_type scalar, const TinyVector& X)noexcept
    {
        return TinyVector(scalar * X.x);
    }

    //! y = x * scalar
    friend TinyVector operator * (const TinyVector& X, value_type scalar)noexcept
    {
        return TinyVector(scalar * X.x);
    }

    //! y = x / scalar
    friend TinyVector operator / (const TinyVector& X, value_type scalar)noexcept
    {
        return TinyVector(X.x / scalar);
    }

    //! y = -x
    friend inline TinyVector operator - (const TinyVector& x)noexcept
    {
        return TinyVector(-x.x);        
    }

    //----------------------------------------
    // other friend method
    //----------------------------------------

    friend value_type dot(const TinyVector& a, const TinyVector& b)noexcept { return a.x * b.x; }

    friend value_type distance   (const TinyVector& a, const TinyVector& b)noexcept { return std::abs(b.x - a.x); }

    friend value_type distance_sq(const TinyVector& a, const TinyVector& b)noexcept { return (b.x - a.x) * (b.x - a.x); }

    friend value_type norm(const TinyVector& a)noexcept { return std::abs(a.x); }

    friend value_type norm_sq(const TinyVector& a)noexcept { return a.x * a.x; }

    friend value_type normalize(TinyVector& a)noexcept
    {
        value_type l = a.x;
        a.x = std::copysign(value_type{ 1 }, a.x);
        return l;
    }

    friend void swap(TinyVector& a, TinyVector& b)noexcept { std::swap(a.x, b.x); }
};

template<typename T>
struct TinyVector<T, 2>
{
    using value_type             = T;
    using size_type              = size_t;
    using reference              = T&;
    using const_reference        = const std::remove_const_t<T>&;
    using pointer                = T*;
    using const_pointer          = const std::remove_const_t<T>*;
    using iterator               = typename StdArray<T, 2>::iterator;
    using const_iterator         = typename StdArray<T, 2>::const_iterator;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    
    value_type x{ 0 }, y{ 0 };

    //----------------------------------------
    // constructor
    //----------------------------------------

    TinyVector() = default;
    TinyVector(const TinyVector&) = default;
    TinyVector& operator = (const TinyVector&) = default;

    TinyVector(value_type X, value_type Y)noexcept :x(X), y(Y) {}
    TinyVector(const value_type (&a)[2])noexcept :x(a[0]), y(a[1]) {}
    TinyVector(const StdArray<T,2>& a)noexcept :x(a[0]), y(a[1]) {}

    //----------------------------------------
    // viewer functions
    //----------------------------------------

    //! @brief create TinyVector object from data pointer.
    static       TinyVector& view(pointer data, size_t offset = 0)noexcept
    {
        //static_assert(sizeof(TinyVector) == sizeof(value_type) * 2, "TinyVector::view(), type size not agree.");
        return *reinterpret_cast<TinyVector*>(data + offset);
    }

    //! @brief create TinyVector object from data pointer.
    static const TinyVector& view(const_pointer data, size_t offset = 0)noexcept
    {
        //static_assert(sizeof(TinyVector) == sizeof(value_type) * 2, "TinyVector::view(), type size not agree.");
        return *reinterpret_cast<const TinyVector*>(data + offset);
    }

    //----------------------------------------
    // common methods
    //----------------------------------------

    constexpr size_type size    ()const noexcept { return 2; }
    constexpr size_type max_size() const noexcept { return 2; }
    constexpr bool      empty   () const noexcept { return false; }

    pointer       data()      noexcept { return &x; }
    const_pointer data()const noexcept { return &x; }

    void assign(value_type X, value_type Y)noexcept { x = X; y = Y; }

    void fill(value_type val)noexcept { x = y = val; }

    void swap(TinyVector& b)noexcept { std::swap(x, b.x); std::swap(y, b.y); }

    reference       max_element()noexcept
    {
        return x >= y ? x : y;
    }
    const_reference max_element()const noexcept
    {
        return x >= y ? x : y;
    }
    reference       min_element()noexcept
    {
        return x <= y ? x : y;
    }
    const_reference min_element()const noexcept
    {
        return x <= y ? x : y;
    }

    //----------------------------------------
    // mathematic methods
    //----------------------------------------

    value_type norm()const noexcept
    {
        return std::hypot(x, y);
    }

    value_type norm_sq()const noexcept
    {
        return x * x + y * y;
    }

    value_type normalize()noexcept
    {
        value_type l = this->norm_sq();
        x /= l != 0 ? l : 1;
        y /= l != 0 ? l : 1;
        return l;
    }

    value_type dot(const TinyVector& b)const noexcept
    {
        return x * b.x + y * b.y;
    }

    value_type cross(const TinyVector& b)const noexcept
    {
        return x * b.y - y * b.x;
    }

    //! Y = a * X + Y
    TinyVector& add_ax(value_type a, const TinyVector& X)noexcept
    {
        x += a * X.x;
        y += a * X.y;
        return *this;
    }

    value_type distance_sq(const TinyVector& b)const noexcept
    {
        auto dx = x - b.x;
        auto dy = y - b.y;
        return dx * dx + dy * dy;
    }

    value_type distance(const TinyVector& b)const noexcept
    {
        auto dx = x - b.x;
        auto dy = y - b.y;
        return std::hypot(dx, dy);
    }

    //----------------------------------------
    // subscript operator
    //----------------------------------------

    reference       operator[](size_type i)      noexcept { ASSERT(i >= 0 && i < size()); return (&x)[i]; }
    const_reference operator[](size_type i)const noexcept { ASSERT(i >= 0 && i < size()); return (&x)[i]; }

    reference       at(size_type i)
    {
        if (i < 0 || i >= size())throw std::out_of_range("TinyVector<*,2>::at(), index out of range");
        return (&x)[i];
    }
    const_reference at(size_type i)const
    {
        if (i < 0 || i >= size())throw std::out_of_range("TinyVector<*,2>::at(), index out of range");
        return (&x)[i];
    }

    template<size_type I>       value_type& entity()      noexcept { static_assert(I >= 0 && I < size(), "TinyVector<*,2>::entity(), index out of range"); return (&x)[I]; }
    template<size_type I> const value_type& entity()const noexcept { static_assert(I >= 0 && I < size(), "TinyVector<*,2>::entity(), index out of range"); return (&x)[I]; }

    //----------------------------------------
    // iterator
    //----------------------------------------

    iterator               begin  ()       noexcept { return (iterator(&x, 0)); }
    const_iterator         begin  () const noexcept { return (const_iterator(&x, 0)); }
    iterator               end    ()       noexcept { return (iterator(&x, size())); }
    const_iterator         end    () const noexcept { return (const_iterator(&x, size())); }
    reverse_iterator       rbegin ()       noexcept { return (reverse_iterator(end())); }
    const_reverse_iterator rbegin () const noexcept { return (const_reverse_iterator(end())); }
    reverse_iterator       rend   ()       noexcept { return (reverse_iterator(begin())); }
    const_reverse_iterator rend   () const noexcept { return (const_reverse_iterator(begin())); }
    const_iterator         cbegin () const noexcept { return (begin()); }
    const_iterator         cend   () const noexcept { return (end()); }
    const_reverse_iterator crbegin() const noexcept { return (rbegin()); }
    const_reverse_iterator crend  () const noexcept { return (rend()); }

    //----------------------------------------
    // assignment operators
    //----------------------------------------
    
    //! assign operator by value list
    template<typename T2>
    TinyVector& operator = (std::initializer_list<T2> list)
    {
        if (list.size() != size())throw std::out_of_range("TinyVector<*,2>::operator=(list), length not agree");
        x = static_cast<value_type>(*list.begin());
        y = static_cast<value_type>(*(list.begin()+1));
        return *this;
    }

    //! assign operator by c array
    template<typename T2>
    TinyVector& operator = (const T2(&v)[2])noexcept
    {
        x = static_cast<value_type>(v[0]);
        y = static_cast<value_type>(v[1]);
        return *this;
    }

    //! assign operator by std::array
    TinyVector& operator = (const std::array<value_type, 2>& v)noexcept
    {
        x = v[0];
        y = v[1];
        return *this;
    }

    //! assign operator by value
    //inline TinyVector& operator = (const value_type& val)
    //{
    //    x = y = val;
    //    return *this;
    //}

    //----------------------------------------
    // mathematic operators
    //----------------------------------------

    bool operator ==(const TinyVector& b)const noexcept
    {
        return x == b.x && y = b.y;
    }

    bool operator !=(const TinyVector& b)const noexcept
    {
        return x != b.x || y != b.y;;
    }

    //! x = scalar * x
    TinyVector& operator *= (value_type scale)noexcept
    {
        x *= scale;
        y *= scale;
        return *this;
    }

    //! x = x / scalar
    TinyVector& operator /= (value_type divisor)noexcept
    {
        x /= divisor;
        y /= divisor;
        return *this;
    }

    //! x = x + y
    TinyVector& operator += (const TinyVector& Y)noexcept
    {
        x += Y.x;
        y += Y.y;
        return *this;
    }

    //! x = x - y
    TinyVector& operator -= (const TinyVector& Y)noexcept
    {
        x -= Y.x;
        y -= Y.y;
        return *this;
    }

    //! z = x + y
    friend TinyVector operator + (const TinyVector& X, const TinyVector& Y)noexcept
    {
        return TinyVector(X.x + Y.x, X.y + Y.y);
    }

    //! z = x - y
    friend TinyVector operator - (const TinyVector& X, const TinyVector& Y)noexcept
    {
        return TinyVector(X.x - Y.x, X.y - Y.y);
    }

    //! y = scalar * x
    friend TinyVector operator * (value_type scalar, const TinyVector& X)noexcept
    {
        return TinyVector(scalar * X.x, scalar * X.y);
    }

    //! y = x * scalar
    friend TinyVector operator * (const TinyVector& X, value_type scalar)noexcept
    {
        return TinyVector(scalar * X.x, scalar * X.y);
    }

    //! y = x / scalar
    friend TinyVector operator / (const TinyVector& X, value_type scalar)noexcept
    {
        return TinyVector(X.x / scalar, X.y / scalar);
    }

    //! y = -x
    friend inline TinyVector operator - (const TinyVector& x)noexcept
    {
        return TinyVector(-x.x, -x.y);
    }
    
    //----------------------------------------
    // friend mathematic functions
    //----------------------------------------

    friend value_type dot(const TinyVector& a, const TinyVector& b)noexcept
    {
        return a.x * b.x + a.y * b.y;
    }

    friend value_type cross(const TinyVector& a, const TinyVector& b)noexcept
    {
        return a.x * b.y - a.y * b.x;
    }

    friend value_type distance(const TinyVector& a, const TinyVector& b)noexcept
    {
        auto dx = b.x - a.x;
        auto dy = b.y - a.y;
        return std::hypot(dx, dy);
    }

    friend value_type distance_sq(const TinyVector& a, const TinyVector& b)noexcept
    {
        auto dx = b.x - a.x;
        auto dy = b.y - a.y;
        return dx * dx + dy * dy;
    }

    friend value_type norm(const TinyVector& a)noexcept
    {
        return std::hypot(a.x, a.y);
    }

    friend value_type norm_sq(const TinyVector& a)noexcept
    {
        return a.x * a.x + a.y * a.y;
    }

    friend value_type normalize(TinyVector& a)noexcept
    {
        value_type l = std::hypot(a.x, a.y);
        a.x /= l != 0 ? l : 1;
        a.y /= l != 0 ? l : 1;
        return l;
    }
};

template<typename T>
struct TinyVector<T, 3>
{
    using value_type             = T;
    using size_type              = size_t;
    using reference              = T&;
    using const_reference        = const std::remove_const_t<T>&;
    using pointer                = T*;
    using const_pointer          = const std::remove_const_t<T>*;
    using iterator               = typename StdArray<T, 3>::iterator;
    using const_iterator         = typename StdArray<T, 3>::const_iterator;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    
    value_type x{ 0 }, y{ 0 }, z{ 0 };

    //----------------------------------------
    // constructor
    //----------------------------------------

    TinyVector() = default;
    TinyVector(const TinyVector&) = default;
    TinyVector& operator = (const TinyVector&) = default;

    TinyVector(value_type X, value_type Y, value_type Z)noexcept :x(X), y(Y), z(Z) {}
    TinyVector(const value_type(&a)[3])noexcept :x(a[0]), y(a[1]), z(a[2]) {}
    TinyVector(const StdArray<T, 3>& a)noexcept :x(a[0]), y(a[1]), z(a[2]) {}

    //----------------------------------------
    // viewer functions
    //----------------------------------------

    //! @brief create TinyVector object from data pointer.
    static TinyVector& view(value_type* data, size_type offset = 0)noexcept
    {
        //static_assert(sizeof(TinyVector) == sizeof(value_type) * 3, "TinyVector::view(), type size not agree.");
        return *reinterpret_cast<TinyVector*>(data + offset);
    }

    //! @brief create TinyVector object from data pointer.
    static const TinyVector& view(const value_type* data, size_type offset = 0)noexcept
    {
        //static_assert(sizeof(TinyVector) == sizeof(value_type) * 3, "TinyVector::view(), type size not agree.");
        return *reinterpret_cast<const TinyVector*>(data + offset);
    }

    //----------------------------------------
    // common methods
    //----------------------------------------

    void assign(value_type X, value_type Y, value_type Z)noexcept { x = X; y = Y; z = Z; }

    constexpr size_type size    ()const noexcept { return 3; }
    constexpr size_type max_size() const noexcept { return 3; }
    constexpr bool      empty   () const noexcept { return false; }

    pointer       data()      noexcept { return &x; }
    const_pointer data()const noexcept { return &x; }

    void fill(value_type val)noexcept { x = y = z = val; }

    void swap(TinyVector& b)noexcept { std::swap(x, b.x); std::swap(y, b.y); std::swap(z, b.z); }

    reference       max_element()noexcept
    {
        return x >= y ? (x >= z ? x : z) : (y >= z ? y : z);
    }
    const_reference max_element()const noexcept
    {
        return x >= y ? (x >= z ? x : z) : (y >= z ? y : z);
    }
    reference       min_element()noexcept
    {
        return x <= y ? (x <= z ? x : z) : (y <= z ? y : z);
    }
    const_reference min_element()const noexcept
    {
        return x <= y ? (x <= z ? x : z) : (y <= z ? y : z);
    }

    //----------------------------------------
    // mathematic methods
    //----------------------------------------

    value_type norm()const noexcept
    {
        return std::hypot(x, y, z);
    }

    value_type norm_sq()const noexcept
    {
        return x * x + y * y + z * z;
    }

    value_type normalize() noexcept
    {
        value_type l = this->norm_sq();
        x /= l != 0 ? l : 1;
        y /= l != 0 ? l : 1;
        z /= l != 0 ? l : 1;
        return l;
    }

    value_type dot(const TinyVector& b)const noexcept
    {
        return x * b.x + y * b.y + z * b.z;
    }

    TinyVector cross(const TinyVector& b)const noexcept
    {
        return TinyVector(
            y * b.z - z * b.y,
            z * b.x - x * b.z,
            x * b.y - y * b.x
        );
    }

    //! @brief Y = a * X + Y 
    TinyVector& add_ax(value_type a, const TinyVector& X)noexcept
    {
        x += a * X.x;
        y += a * X.y;
        z += a * X.z;
        return *this;
    }

    value_type distance_sq(const TinyVector& b)const noexcept
    {
        auto dx = x - b.x;
        auto dy = y - b.y;
        auto dz = z - b.z;
        return dx * dx + dy * dy + dz * dz;
    }

    value_type distance(const TinyVector& b)const noexcept
    {
        auto dx = x - b.x;
        auto dy = y - b.y;
        auto dz = z - b.z;
        return std::hypot(dx, dy, dz);
    }

    //----------------------------------------
    // subscript operator
    //----------------------------------------

    reference       operator[](size_type i)      noexcept { ASSERT(i >= 0 && i < size()); return (&x)[i]; }
    const_reference operator[](size_type i)const noexcept { ASSERT(i >= 0 && i < size()); return (&x)[i]; }

    reference       at(size_type i)
    {
        if (i < 0 || i >= size())throw std::out_of_range("TinyVector<*,3>::at(), index out of range");
        return (&x)[i];
    }
    const_reference at(size_type i)const
    {
        if (i < 0 || i >= size())throw std::out_of_range("TinyVector<*,3>::at(), index out of range");
        return (&x)[i];
    }

    template<size_type I>       value_type& entity()      noexcept { static_assert(I < size(), "TinyVector<*,3>::entity(), index out of range"); return (&x)[I]; }
    template<size_type I> const value_type& entity()const noexcept { static_assert(I < size(), "TinyVector<*,3>::entity(), index out of range"); return (&x)[I]; }

    //----------------------------------------
    // iterator
    //----------------------------------------
                                  
    iterator               begin  ()       noexcept { return (iterator(&x, 0)); }
    const_iterator         begin  () const noexcept { return (const_iterator(&x, 0)); }
    iterator               end    ()       noexcept { return (iterator(&x, size())); }
    const_iterator         end    () const noexcept { return (const_iterator(&x, size())); }
    reverse_iterator       rbegin ()       noexcept { return (reverse_iterator(end())); }
    const_reverse_iterator rbegin () const noexcept { return (const_reverse_iterator(end())); }
    reverse_iterator       rend   ()       noexcept { return (reverse_iterator(begin())); }
    const_reverse_iterator rend   () const noexcept { return (const_reverse_iterator(begin())); }
    const_iterator         cbegin () const noexcept { return (begin()); }
    const_iterator         cend   () const noexcept { return (end()); }
    const_reverse_iterator crbegin() const noexcept { return (rbegin()); }
    const_reverse_iterator crend  () const noexcept { return (rend()); }

    //----------------------------------------
    // assignment operator
    //----------------------------------------
    
    //!assign operator by value list
    template<typename T2>
    TinyVector& operator = (std::initializer_list<T2> list)
    {
        if (list.size() != size())throw std::out_of_range("TinyVector<*,3>::operator=(list), length not agree");
        x = static_cast<value_type>(*list.begin());
        y = static_cast<value_type>(*(list.begin() + 1));
        z = static_cast<value_type>(*(list.begin() + 2));
        return *this;
    }

    //!assign operator by c array
    template<typename T2>
    TinyVector& operator = (const T2(&v)[3])noexcept
    {
        x = static_cast<value_type>(v[0]);
        y = static_cast<value_type>(v[1]);
        z = static_cast<value_type>(v[2]);
        return *this;
    }

    //! assign operator by std::array
    TinyVector& operator = (const std::array<value_type, 3>& v)noexcept
    {
        x = v[0];
        y = v[1];
        z = v[2];
        return *this;
    }

    //! assign operator by value
    //TinyVector& operator = (const value_type& val)
    //{
    //    x = y = z = val;
    //    return *this;
    //}

    //----------------------------------------
    // mathematic operator
    //----------------------------------------

    bool operator ==(const TinyVector& b)const noexcept
    {
        return x == b.x && y == b.y && z == b.z;
    }

    bool operator !=(const TinyVector& b)const noexcept
    {
        return x != b.x || y != b.y || z != b.z;
    }

    //! x = scalar * x
    TinyVector& operator *= (const value_type& scale)noexcept
    {
        x *= scale;
        y *= scale;
        z *= scale;
        return *this;
    }

    //! x = x / scalar
    TinyVector& operator /= (const value_type& divisor)noexcept
    {
        x /= divisor;
        y /= divisor;
        z /= divisor;
        return *this;
    }

    //! x = x + y
    TinyVector& operator += (const TinyVector& Y)noexcept
    {
        x += Y.x;
        y += Y.y;
        z += Y.z;
        return *this;
    }

    //! x = x - y
    TinyVector& operator -= (const TinyVector& Y)noexcept
    {
        x -= Y.x;
        y -= Y.y;
        z -= Y.z;
        return *this;
    }

    //! z = x + y
    friend TinyVector operator + (const TinyVector& X, const TinyVector& Y)noexcept
    {
        return TinyVector(X.x + Y.x, X.y + Y.y, X.z + Y.z);
    }

    //! z = x - y
    friend TinyVector operator - (const TinyVector& X, const TinyVector& Y)noexcept
    {
        return TinyVector(X.x - Y.x, X.y - Y.y, X.z - Y.z);
    }

    //! y = scalar * x
    friend TinyVector operator * (value_type scalar, const TinyVector& X)noexcept
    {
        return TinyVector(scalar * X.x, scalar * X.y, scalar * X.z);
    }

    //! y = x * scalar
    friend TinyVector operator * (const TinyVector& X, value_type scalar)noexcept
    {
        return TinyVector(scalar * X.x, scalar * X.y, scalar * X.z);
    }

    //! y = x / scalar
    friend TinyVector operator / (const TinyVector& X, value_type scalar)noexcept
    {
        return TinyVector(X.x / scalar, X.y / scalar, X.z / scalar);
    }

    //! y = -x
    friend TinyVector operator - (const TinyVector& x)noexcept
    {
        return TinyVector(-x.x, -x.y, -x.z);
    }

    //----------------------------------------
    // friend mathematic functions
    //----------------------------------------

    friend value_type dot(const TinyVector& a, const TinyVector& b)noexcept
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    friend TinyVector cross(const TinyVector& a, const TinyVector& b)noexcept
    {
        return TinyVector(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        );
    }

    friend value_type distance(const TinyVector& a, const TinyVector& b)noexcept
    {
        auto dx = b.x - a.x;
        auto dy = b.y - a.y;
        auto dz = b.z - a.z;
        return std::hypot(dx, dy, dz);
    }

    friend value_type distance_sq(const TinyVector& a, const TinyVector& b)noexcept
    {
        auto dx = b.x - a.x;
        auto dy = b.y - a.y;
        auto dz = b.z - a.z;
        return dx * dx + dy * dy + dz * dz;
    }

    friend value_type norm(const TinyVector& a)noexcept
    {
        return std::hypot(a.x, a.y, a.z);
    }

    friend value_type norm_sq(const TinyVector& a)noexcept
    {
        return a.x * a.x + a.y * a.y + a.z * a.z;
    }

    friend value_type normalize(TinyVector& a)noexcept
    {
        value_type l = std::hypot(a.x, a.y, a.z);
        a.x /= l != 0 ? l : 1;
        a.y /= l != 0 ? l : 1;
        a.z /= l != 0 ? l : 1;
        return l;
    }
};

using Vec2 = TinyVector<double, 2>;
using Vec3 = TinyVector<double, 3>;
