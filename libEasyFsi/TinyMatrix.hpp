#pragma once
#include <array>
#include <initializer_list>

#include "TinyVector.hpp"

//! row-major fixed-size matrix class
//! @param T  the element value type.
//! @param M  the row count.
//! @param N  the column count, default is M.
template<typename T, size_t M, size_t N = M>
class TinyMatrix
{
public:
    using value_type             = T;     //!< the element value type
    using size_type              = size_t;//!< the size/indice value type
    using pointer                = T* ;
    using const_pointer          = const std::remove_const_t<T>*;
    using reference              = T & ;
    using const_reference        = const std::remove_const_t<T>&;
    using difference_type        = ptrdiff_t;
    using iterator               = typename std::array<T, M * N>::iterator;
    using const_iterator         = typename std::array<T, M * N>::const_iterator;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    static_assert(M > 0, "TinyMatrix<T, M, N>, matrix row count must great than zero"   );
    static_assert(N > 0, "TinyMatrix<T, M, N>, matrix column count must great than zero");

    static constexpr size_t dim_min = M <= N ? M : N;//!<the minimum dimension

    TinyMatrix() = default;
    TinyMatrix(const TinyMatrix&) = default;
    TinyMatrix& operator = (const TinyMatrix&) = default;

    explicit TinyMatrix(std::initializer_list<std::initializer_list<value_type>> list)
    {
        static_assert(M > 0 && N > 0, "TinyMatrix(), matrix dimensions must great than zero");

        auto it = list.begin();
        for (size_t i = 0; i < M && i < list.size(); ++i, ++it) {
            auto it2 = (*it).begin();
            for (size_t j = 0; j < N && it2 != (*it).end(); ++j, ++it2) {
                data_[i * N + j] = *it2;
            }
        }
    }

    TinyMatrix& operator = (std::initializer_list<std::initializer_list<value_type> > list)
    {
        auto it = list.begin();
        for (size_t i = 0; i < M && i < list.size(); ++i, ++it) {
            auto it2 = (*it).begin();
            for (size_t j = 0; j < N && it2 != (*it).end(); ++j, ++it2) {
                data_[i * N + j] = *it2;
            }
        }
        return *this;
    }

    constexpr size_t row_count   () const noexcept { return M; }
    constexpr size_t column_count() const noexcept { return N; }
    constexpr size_type max_size() const noexcept { return M * N; }
    constexpr bool is_squared() const noexcept { return M == N; }
    constexpr bool empty() const noexcept { return false; }
    constexpr bool row_major()const { return true; }
    constexpr bool column_major()const { return false; }
    constexpr bool contiguous()const { return true; }
    constexpr size_type numel() const noexcept { return M * N; }

    pointer       data()      noexcept { return data_; }
    const_pointer data()const noexcept { return data_; }

    size_type size(size_t i)const { ASSERT(i >= 0 && i < 2); return i = 0 ? M : N; }

    std::array<size_t, 2> size()const noexcept { return std::array<size_t, 2>{M, N}; }

    //! get element value by 2d subscript, \i = row index, \j = column index
    reference operator()(size_t i, size_t j)
    {
        ASSERT(i >= 0 && i < M && j >= 0 && j < N);
        return data_[i * N + j];
    }
    //! get element value by 2d subscript, \i = row index, \j = column index
    const_reference operator()(size_t i, size_t j)const
    {
        ASSERT(i >= 0 && i < M && j >= 0 && j < N);
        return data_[i * N + j];
    }

    //! get element value by 1d subscript (row-major)
    reference operator()(size_t i)
    {
        ASSERT(i >= 0 && i < M * N);
        return data_[i];
    }
    //! get element value by 1d subscript (row-major)
    const_reference operator()(size_t i)const
    {
        ASSERT(i >= 0 && i < M * N);
        return data_[i];
    }

    //! get i-row of matrix, return reference a vector
    //TinyVector<T, N>& operator[](size_t i)
    //{
    //    ASSERT(i >= 0 && i < M);
    //    return TinyVector<T, N>::view(data_ + i * N);
    //}
    ////!get i-row of matrix, return reference a vector
    //const TinyVector<T, N>& operator[](size_t i)const
    //{
    //    ASSERT(i >= 0 && i < M);
    //    return TinyVector<T, N>::view(data_ + i * N);
    //}

    //! get one row using zero-based index, is same with opreator []
    TinyVector<T, N>& row(size_t i)
    {
        ASSERT(i >= 0 && i < M);
        return TinyVector<T, N>::view(data_ + i * N);
    }
    //! get one row using zero-based index, is same with opreator []
    const TinyVector<T, N>& row(size_t i)const
    {
        ASSERT(i >= 0 && i < M);
        return TinyVector<T, N>::view(data_ + i * N);
    }

    //! get element value by 2d subscript
    template<size_t I, size_t J>       value_type& entity()
    {
        static_assert(I < M && J < N, "TinyMatrix::entity<I,J>(), subscript out of range.");
        return data_[I * N + J];
    }
    //! get element value by 2d subscript
    template<size_t I, size_t J> const value_type& entity()const
    {
        static_assert(I < M && J < N, "TinyMatrix::entity<I,J>(), subscript out of range.");
        return data_[I * N + J];
    }

    //-------------------------------------------
    // iterator
    //-------------------------------------------

    iterator               begin  ()      noexcept { return iterator(data_); }
    iterator               end    ()      noexcept { return iterator(data_, M * N); }
    const_iterator         begin  ()const noexcept { return const_iterator(data_); }
    const_iterator         end    ()const noexcept { return const_iterator(data_, M * N); }
    const_iterator         cbegin ()const noexcept { return const_iterator(data_); }
    const_iterator         cend   ()const noexcept { return const_iterator(data_, M * N); }
    reverse_iterator       rbegin ()      noexcept { return (reverse_iterator(end())); }
    const_reverse_iterator rbegin ()const noexcept { return (const_reverse_iterator(end())); }
    reverse_iterator       rend   ()      noexcept { return (reverse_iterator(begin())); }
    const_reverse_iterator rend   ()const noexcept { return (const_reverse_iterator(begin())); }
    const_reverse_iterator crbegin()const noexcept { return (rbegin()); }
    const_reverse_iterator crend  ()const noexcept { return (rend()); }

    //! @brief get element using zero-based index
    //! @param i  the zero-based element index
    reference at(size_t i)
    {
        if (i < 0 || i >= M * N)
            throw std::out_of_range("TinyMatrix::at(), index out of range.");
        return data_[i];
    }
    //! @brief get element using zero-based index
    //! @param i  the zero-based element index
    const_reference at(size_t i)const
    {
        if (i < 0 || i >= M * N)
            throw std::out_of_range("TinyMatrix::at(), index out of range.");
        return data_[i];
    }

    reference       front()      noexcept { return *data_; }
    const_reference front()const noexcept { return *data_; }
    reference       back ()      noexcept { return *(data_ + M * N - 1); }
    const_reference back ()const noexcept { return *(data_ + M * N - 1); }

    //-------------------------------------------
    // arithmetic operator
    //-------------------------------------------

    //! A = A * scalar
    TinyMatrix& operator *= (const value_type& scale)noexcept
    {
        for (size_t i = 0; i < M * N; ++i) data_[i] *= scale;
        return *this;
    }
    //! A = A / scalar
    TinyMatrix& operator /= (const value_type& divisor)noexcept
    {
        for (size_t i = 0; i < M * N; ++i) data_[i] /= divisor;
        return *this;
    }
    //! A = A + m
    TinyMatrix& operator += (const TinyMatrix& m)noexcept
    {
        for (size_t i = 0; i < M * N; ++i) data_[i] += m.data_[i];
        return *this;
    }
    //! A = A - m
    TinyMatrix& operator -= (const TinyMatrix& m)noexcept
    {
        for (size_t i = 0; i < M * N; ++i) data_[i] -= m.data_[i];
        return *this;
    }
    //! A = A * m
    TinyMatrix& operator *= (const TinyMatrix<T, N, N>& m)noexcept
    {
        TinyVector<T, N> ri;//the i-row

        //loop each row
        for (size_t i = 0; i < M; ++i) {
            //copy i-row to rj
            auto ai = data_ + i * N;
            for (size_t j = 0; j < N; ++j)ri[j] = ai[j];

            //loop each column
            for (size_t j = 0; j < N; ++j) {
                ai[j] = 0;
                // Aik * mkj
                for (size_t k = 0; k < N; ++k) {
                    ai[j] += ri[k] * m(k, j);
                }
            }
        }
        return *this;
    }

    //! @brief 将所有元素置为0
    TinyMatrix& zero()noexcept
    {
        for (size_t i = 0; i < M * N; ++i) data_[i] = 0;
        return *this;
    }

    //! @brief 将矩阵对角线元素置为1，非对角线置为0
    TinyMatrix& identity()noexcept
    {
        zero();
        for (size_t i = 0; i < dim_min; ++i)data_[i * N + i] = 1;
        return *this;
    }

    //! @brief 将矩阵所有元素置为给定值
    TinyMatrix& fill(const value_type& value)noexcept
    {
        for (size_t i = 0; i < M * N; ++i) data_[i] = value;
        return *this;
    }

    //! @brief 对矩阵进行转置操作
    TinyMatrix& transpose()noexcept
    {
        static_assert(M == N, "matrix is not squared");

        for (size_t i = 0; i < M; ++i) {
            for (size_t j = i + 1; j < M; ++j) {
                std::swap(data_[i * M + j], data_[j * M + i]);
            }
        }
        return *this;
    }

    //! @brief 设置为指定矩阵的转置
    template<typename T2 = T>
    TinyMatrix& transpose_from(const TinyMatrix<T2, N, M>& m)noexcept
    {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                data_[i * N + j] = m(j, i);
            }
        }
        return *this;
    }

    //! @brief 设置为指定矩阵的转置
    template<typename T2 = T>
    TinyMatrix& transpose_from(const T2(&m)[N][M])noexcept
    {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                data_[i * N + j] = m[j][i];
            }
        }
        return *this;
    }

    //! @brief 设置为指定矩阵的转置
    template<typename MAT_NM>
    TinyMatrix& transpose_from(const MAT_NM& m)
    {
        if (m.row_count() != N || m.column_count() != M)
            throw std::out_of_range("TinyMatrix::transpose_from(), dimension not agree.");

        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                data_[i * N + j] = m(j, i);
            }
        }
        return *this;
    }

    //! A = A * m
    TinyMatrix& right_multiple(const TinyMatrix<T, N, N>& m)noexcept
    {
        return this->operator*=(m);
    }

    //! A = m * A
    TinyMatrix& left_multiply(const TinyMatrix<T, M, M>& m)noexcept
    {
        TinyVector<T, M> cj;//the j-column

        //loop each column
        for (size_t j = 0; j < N; ++j) {
            //copy j-column to cj
            for (size_t i = 0; i < M; ++i)cj[i] = data_[i * N + j];

            //loop each row
            for (size_t i = 0; i < M; ++i) {
                data_[i * N + j] = dot(m.row(i), cj);
            }
        }
        return *this;
    }

    // y = A * x
    void apply(const TinyVector<T, N>& x, TinyVector<T, M>& y)const noexcept
    {
        for (size_t i = 0; i < M; ++i) {
            auto ri = data_ + i * N;//the i-row
            y[i] = 0;
            for (size_t j = 0; j < N; ++j)y[i] += ri[j] * x[j];
        }
    }
    // v = A * v
    template<typename VEC>
    void apply(VEC& v)const
    {
        static_assert(M == N, "TinyMatrix::ApplyTo(vec), matrix is not square.");
        TinyVector<T, N> x;
        for (size_t i = 0; i < N; ++i)x[i] = v[i];
        apply(x, v);
    }

    // C = A * B
    template<size_t L>
    void apply(const TinyMatrix<T, N, L>& B, TinyMatrix<T, M, L>& C)const noexcept
    {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < L; ++j) {
                C(i, j) = 0;
                for (size_t k = 0; k < N; ++k) {
                    C(i, j) += data_[i * N + k] * B(k, j);
                }
            }
        }
    }
    
    // C = C + A * B
    template<size_t L>
    void apply_add(const TinyMatrix<T, N, L>& B, TinyMatrix<T, M, L>& C)const noexcept
    {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < L; ++j) {
                for (size_t k = 0; k < N; ++k) {
                    C(i, j) += data_[i * N + k] * B(k, j);
                }
            }
        }
    }
    
    // C = A^T * B
    template<size_t L>
    void transpose_apply(const TinyMatrix<T, M, L>& B, TinyMatrix<T, N, L>& C)const noexcept
    {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < L; ++j) {
                C(i, j) = 0;
                for (size_t k = 0; k < M; ++k) {
                    C(i, j) += data_[k * N + i] * B(k, j);
                }
            }
        }
    }
    
    // C = C + A^T * B
    template<size_t L>
    void transpose_apply_add(const TinyMatrix<T, M, L>& B, TinyMatrix<T, N, L>& C)const noexcept
    {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < L; ++j) {
                for (size_t k = 0; k < M; ++k) {
                    C(i, j) += data_[k * N + i] * B(k, j);
                }
            }
        }
    }
    
    //! @brief create TinyMatrix object from data pointer.
    static TinyMatrix& view(value_type* data, size_t offset = 0)noexcept
    {
        //static_assert(sizeof(TinyMatrix) == sizeof(value_type) * M * N, "TinyMatrix::view_of(), type size not agree.");
        return *reinterpret_cast<TinyMatrix*>(data + offset);
    }

    //! @brief create TinyMatrix object from data pointer.
    static const TinyMatrix& view(const value_type* data, size_t offset = 0)noexcept
    {
        //static_assert(sizeof(TinyMatrix) == sizeof(value_type) * M * N, "TinyMatrix::view_of(), type size not agree.");
        return *reinterpret_cast<const TinyMatrix*>(data + offset);
    }

    //! @brief create TinyMatrix object from data pointer.
    static TinyMatrix& view(value_type (&data)[M][N])noexcept
    {
        //static_assert(sizeof(TinyMatrix) == sizeof(value_type) * M * N, "TinyMatrix::view_of(), type size not agree.");
        return *reinterpret_cast<TinyMatrix*>(&data[0][0]);
    }

    //! @brief create TinyMatrix object from data pointer.
    static const TinyMatrix& view(const value_type (&data)[M][N])noexcept
    {
        //static_assert(sizeof(TinyMatrix) == sizeof(value_type) * M * N, "TinyMatrix::view_of(), type size not agree.");
        return *reinterpret_cast<const TinyMatrix*>(&data[0][0]);
    }

private:
    value_type data_[M * N]{ 0 };
};

template<typename T, size_t M, size_t N>
TinyVector<T, N> dot(const TinyVector<T, M>& a, const TinyMatrix<T, M, N>& b)noexcept
{
    TinyVector<T, N> c;
    for (size_t j = 0; j < N; ++j) {
        c[j] = 0;
        for (size_t i = 0; i < M; ++i)
            c[j] += a[i] * b(i, j);
    }
    return c;
}
template<typename T, size_t M, size_t N>
TinyVector<T, N> dot(const TinyMatrix<T, 1, M>& a, const TinyMatrix<T, M, N>& b)noexcept
{
    TinyVector<T, N> c;
    for (size_t j = 0; j < N; ++j) {
        c[j] = 0;
        for (size_t i = 0; i < M; ++i)
            c[j] += a(i) * b(i, j);
    }
    return c;
}
template<typename T, size_t M, size_t N>
TinyVector<T, M> dot(const TinyMatrix<T, M, N>& a, const TinyVector<T, N>& b)noexcept
{
    TinyVector<T, M> c;
    for (size_t i = 0; i < M; ++i) {
        c[i] = 0;
        for (size_t j = 0; j < N; ++j)
            c[i] += a(i, j) * b[j];
    }
    return c;
}
template<typename T, size_t M, size_t N, size_t L>
TinyMatrix<T, M, L> dot(const TinyMatrix<T, M, N>& a, const TinyMatrix<T, N, L>& b)noexcept
{
    TinyMatrix<T, M, L> c;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < L; ++j) {
            c(i, j) = 0;
            for (size_t k = 0; k < N; ++k)
                c(i, j) += a(i, k) * b(k, j);
        }
    }
    return c;
}

//! @brief 计算 1*1 矩阵的行列式
template<typename T> T det(const TinyMatrix<T, 1, 1>& m)noexcept
{
    return m(0);
}

//! @brief 计算 2*2 矩阵的行列式
template<typename T> T det(const TinyMatrix<T, 2, 2>& m)noexcept
{
    return
          m(0, 0) * m(1, 1)
        - m(0, 1) * m(1, 0);
}

//! @brief 计算 3*3 矩阵的行列式
template<typename T> T det(const TinyMatrix<T, 3, 3>& m)noexcept
{
    return
        m(0, 0) * m(1, 1) * m(2, 2) +
        m(0, 1) * m(1, 2) * m(2, 0) +
        m(0, 2) * m(1, 0) * m(2, 1) -
        m(0, 0) * m(1, 2) * m(2, 1) -
        m(0, 1) * m(1, 0) * m(2, 2) -
        m(0, 2) * m(1, 1) * m(2, 0);
}

//! @brief 计算 4*4 矩阵的行列式
template<typename T> T det(const TinyMatrix<T, 4, 4>& mat)noexcept
{
    auto m = mat.data();
    auto inv0 =
        m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15]
        + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
    auto inv4 =
        -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15]
        - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
    auto inv8 =
        m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15]
        + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
    auto inv12 = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14];
    return m[0] * inv0 + m[1] * inv4 + m[2] * inv8 + m[3] * inv12;
}

//! @brief 计算 1*1 矩阵的逆: m = m^-1
//! @return true=sucessfull, false=singular
template<typename T> bool inverse(TinyMatrix<T, 1, 1>& m)noexcept
{
    if (m(0) != T(0)) {
        m(0) = T(1) / m(0);
        return true;
    }
    else {
        m(0) = 0;
        return false;
    }
}

//! @brief 计算 1*1 矩阵的逆: minv = m^-1
//! @return true=sucessfull, false=singular
template<typename T> bool inverse(const TinyMatrix<T, 1, 1>& m, TinyMatrix<T, 1, 1>& minv)noexcept
{
    if (m(0) != T(0)) {
        minv(0) = T(1) / m(0);
        return true;
    }
    else {
        minv(0) = T(0);
        return false;
    }
}

//! @brief 计算 2*2 矩阵的逆: m = m^-1
//! @return true=sucessfull, false=singular
template<typename T> bool inverse(TinyMatrix<T, 2, 2>& m)noexcept
{
    T d = det(m);

    if (d != T(0)) {
        T m00 = m(0, 0);
        T m01 = m(0, 1);
        m(0, 0) =  m(1, 1) / d;
        m(1, 1) =  m00     / d;
        m(0, 1) = -m(0, 1) / d;
        m(1, 0) = -m01     / d;
        return true;
    }
    else {
        m.zero();
        return false;
    }
}

//! @brief 计算 2*2 矩阵的逆: minv = m^-1
//! @return true=sucessfull, false=singular
template<typename T> bool inverse(const TinyMatrix<T, 2, 2>& m, TinyMatrix<T, 2, 2>& minv)noexcept
{
    if (m.data() == minv.data())return inverse(minv);

    T d = det(m);

    if (d != T(0)) {
        minv(0, 0) =  m(1, 1) / d;
        minv(0, 1) = -m(0, 1) / d;
        minv(1, 0) = -m(1, 0) / d;
        minv(1, 1) =  m(0, 0) / d;
        return true;
    }
    else {
        minv.zero();
        return false;
    }
}

//! @brief 计算 3*3 矩阵的逆: m = m^-1
//! @return true=sucessfull, false=singular
template<typename T> bool inverse(TinyMatrix<T, 3, 3>& m)noexcept
{
    //details::GLmatOp<ROW_MAJOR>::InverseMat3x3(m);
    auto d = m.det();
    if (std::abs(d) == 0) return false;

    T tm[3 * 3];
    static_assert(sizeof(tm) == sizeof(m), "size not agree");
    std::memcpy(tm, m.data(), sizeof(tm));
    auto inv = m.data();
    auto r = T{ 1 } / d;
    inv[0] = r * (tm[4] * tm[8] - tm[5] * tm[7]);
    inv[1] = r * (tm[2] * tm[7] - tm[1] * tm[8]);
    inv[2] = r * (tm[1] * tm[5] - tm[2] * tm[4]);
    inv[3] = r * (tm[5] * tm[6] - tm[3] * tm[8]);
    inv[4] = r * (tm[0] * tm[8] - tm[2] * tm[6]);
    inv[5] = r * (tm[2] * tm[3] - tm[0] * tm[5]);
    inv[6] = r * (tm[3] * tm[7] - tm[4] * tm[6]);
    inv[7] = r * (tm[1] * tm[6] - tm[0] * tm[7]);
    inv[8] = r * (tm[0] * tm[4] - tm[1] * tm[3]);
    return true;
}

//! @brief 计算 3*3 矩阵的逆: minv = m^-1
//! @return true=sucessfull, false=singular
template<typename T> bool inverse(const TinyMatrix<T, 3, 3>& m, TinyMatrix<T, 3, 3>& minv)noexcept
{
    if (m.data() == minv.data())
        return inverse(minv);
    else {
        auto d = m.det();
        if (std::abs(d) == 0) return false;
        auto tm  = m.data();
        auto inv = minv.data();
        auto r = T{ 1 } / d;
        inv[0] = r * (tm[4] * tm[8] - tm[5] * tm[7]);
        inv[1] = r * (tm[2] * tm[7] - tm[1] * tm[8]);
        inv[2] = r * (tm[1] * tm[5] - tm[2] * tm[4]);
        inv[3] = r * (tm[5] * tm[6] - tm[3] * tm[8]);
        inv[4] = r * (tm[0] * tm[8] - tm[2] * tm[6]);
        inv[5] = r * (tm[2] * tm[3] - tm[0] * tm[5]);
        inv[6] = r * (tm[3] * tm[7] - tm[4] * tm[6]);
        inv[7] = r * (tm[1] * tm[6] - tm[0] * tm[7]);
        inv[8] = r * (tm[0] * tm[4] - tm[1] * tm[3]);
        return true;
    }
}

//! @brief 计算 4*4 矩阵的逆
//! @return true=sucessfull, false=singular
template<typename T> bool inverse(TinyMatrix<T, 4, 4>& mat)
{
    auto m = mat.data();
    T tm[4 * 4];

    tm[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15]
        + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
    tm[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15]
        - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
    tm[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15]
        + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
    tm[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14]
        - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
    tm[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15]
        - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
    tm[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15]
        + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
    tm[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15]
        - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
    tm[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14]
        + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
    tm[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15]
        + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
    tm[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15]
        - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
    tm[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15]
        + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
    tm[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14]
        - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
    tm[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11]
        - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
    tm[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11]
        + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
    tm[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11]
        - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
    tm[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10]
        + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

    auto d = m[0] * tm[0] + m[1] * tm[4] + m[2] * tm[8] + m[3] * tm[12];
    auto r = d != 0 ? T{ 1 } / d : T{ 0 };
    for (size_t i = 0; i < 16; i++)m[i] = tm[i] * r;

    return std::abs(d) > 0;
}

//! @brief 计算 4*4 矩阵的逆: minv = m^-1
//! @return true=sucessfull, false=singular
template<typename T> bool inverse(const TinyMatrix<T, 4, 4>& mat, TinyMatrix<T, 4, 4>& minv)noexcept
{
    auto m  = mat.data();
    auto tm = minv.data();
    if (m == tm)return inverse(minv);

    tm[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15]
        + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
    tm[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15]
        - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
    tm[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15]
        + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
    tm[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14]
        - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
    tm[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15]
        - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
    tm[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15]
        + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
    tm[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15]
        - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
    tm[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14]
        + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
    tm[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15]
        + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
    tm[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15]
        - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
    tm[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15]
        + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
    tm[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14]
        - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
    tm[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11]
        - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
    tm[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11]
        + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
    tm[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11]
        - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
    tm[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10]
        + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

    auto d = m[0] * tm[0] + m[1] * tm[4] + m[2] * tm[8] + m[3] * tm[12];
    auto r = d != 0 ? T{ 1 } / d : T{ 0 };
    for (size_t i = 0; i < 16; i++)tm[i] *= r;

    return std::abs(d) > 0;
}

//! @brief 计算 N*N 矩阵的逆: m = m^-1
//! @return true=sucessfull, false=singular
template<typename T, size_t Rank> bool inverse(TinyMatrix<T, Rank, Rank>& m)noexcept
{
    static_assert(Rank > 0, "InverseMatrix(), matrix rank is zero");

    constexpr int rank = static_cast<int>(Rank);

    //---- 全选主元法对矩阵求逆 ----

    int is[rank], js[rank];

    // 循环所有行
    for (int k = 0; k < rank; ++k) {

        // 查找右下矩阵块中的绝对值最大的元素所在行和列
        T d(0);
        for (int i = k; i < rank; ++i) {
            for (int j = k; j < rank; ++j) {
                T p = std::abs(m(i, j));//? must use std::abs here.
                if (p > d) {
                    d = p;
                    is[k] = i;
                    js[k] = j;
                }
            }
        }

        if (d == 0) {
            //throw std::runtime_error("Matrix::inverse() -- singular matrix");
            m.zero();
            return false;
        }

        auto rk = &m(k * rank);//矩阵的第 k 行

        // 如果最大元素所在行不是第 k 行，则与第 k 行交换
        if (is[k] != k) {
            auto rk2 = &m(is[k] * rank);//矩阵的第is[k]行
            for (int j = 0; j < rank; ++j) {
                std::swap(rk[j], rk2[j]);
            }
        }

        // 如果最大元素所在列不是第 k 列，则与第 k 列交换
        if (js[k] != k) {
            int j = js[k];
            for (int i = 0; i < rank; ++i) {
                std::swap(m(i, k), m(i, j));
            }
        }

        // 当前行对角项置为自身倒数
        T rkk = 1 / rk[k];
        // 当前行除以对角项元素，不包括对角项
        for (int j = 0; j < rank; ++j)rk[j] *= rkk;
        rk[k] = rkk;

        //
        for (int i = 0; i<rank; i++) {
            if (i != k) {
                auto ri = &m(i * rank);//矩阵的第i行
                for (int j = 0; j<rank; j++) {
                    if (j != k) {
                        ri[j] -= ri[k] * rk[j];
                    }
                }
            }
        }

        // 第k列乘以 -第k行，不包括对角项
        for (int i = 0; i<rank; i++) {
            if (i != k) {
                m(i, k) *= -rk[k];
            }
        }
    }

    for (int k = rank - 1; k >= 0; --k) {
        T* rk = m.data() + (size_t)k * (size_t)rank;
        if (js[k] != k) {
            T* rk2 = m.data() + (size_t)js[k] * (size_t)rank;
            for (int j = 0; j < rank; ++j) {
                std::swap(rk[j], rk2[j]);
            }
        }
        if (is[k] != k) {
            int k2 = is[k];
            for (int i = 0; i < rank; ++i) {
                std::swap(m(i, k), m(i, k2));
            }
        }
    }

    return true;
}

//! @brief 计算 N*N 矩阵的逆: m = m^-1
//! @return true=sucessfull, false=singular
template<typename T, size_t Rank> bool inverse(const TinyMatrix<T, Rank, Rank>& M, TinyMatrix<T, Rank, Rank>& invM)noexcept
{
    if (M.data() != invM.data())invM = M; // copy firstly
    return inverse(invM);    
}

//! create a unit matrix, is same with MATLAB eye function
template<typename T, size_t M, size_t N = M> TinyMatrix<T, M, N> eye()noexcept
{
    return TinyMatrix<T, M, N>().identity();
}

//! create a matrix, is same with MATLAB zeros function
template<typename T, size_t M, size_t N = M> TinyMatrix<T, M, N> zeros()noexcept
{
    return TinyMatrix<T, M, N>().zeros();
}

//! create a matrix, is same with MATLAB ones function
template<typename T, size_t M, size_t N = M> TinyMatrix<T, M, N> ones()noexcept
{
    return TinyMatrix<T, M, N>().fill(1);
}

//! create a diagonal matrix by given vector, is same with MATLAB diag function
template<typename T, size_t N> TinyMatrix<T, N, N> diag(const TinyVector<T, N>& vec)noexcept
{
    TinyMatrix<T, N, N> m;
    m.zero();
    for (size_t i = 0; i < N; ++i)m(i, i) = vec[i];
    return m;
}

//! generate diagonal vector of matrix, is same with MATLAB diag function
template<typename T, size_t N> TinyVector<T, N> diag(const TinyMatrix<T, N, N>& m)noexcept
{
    TinyVector<T, N> v;
    for (size_t i = 0; i < N; ++i)v[i] = m(i, i);
    return v;
}

//! generate transposed matrix
template<typename T, size_t M, size_t N = M> TinyMatrix<T, N, M> transpose(const TinyMatrix<T, M, N>& m)noexcept
{
    TinyMatrix<T, N, M> mt;
    mt.transpose_from(m);
    return mt;
}

//! generate transposed matrix
template<typename T, size_t M, size_t N = M> void transpose(const TinyMatrix<T, M, N>& m, TinyMatrix<T, N, M>& mt)noexcept
{
    mt.transpose_from(m);
    return mt;
}

template<typename T, size_t M>
TinyVector<T, M> transpose(const TinyMatrix<T, M, 1>& a)noexcept
{
    TinyMatrix<T, M> b;
    for (size_t i = 0; i < M; ++i)
        b(i) = a(i);
    return b;
}

/*
void Test(bool condition, const char* msg)
{
    if (condition) {
        printf("[PASS] %s\n", msg);
    }
    else {
        printf("[FAIL] %s\n", msg);
    }
}
void TestTinyMatrix()
{
    double a[2] = { 1, 2 };
    double b[2][2] = { {1, 2}, {3, 4} };

    auto& x = TinyVector<double, 2>::view_of(a);
    TinyVector<double, 2> y;

    auto& A = TinyMatrix<double, 2, 2>::view_of(&b[0][0]);
    Test(A(0, 0) == 1 && A(0, 1) == 2 && A(1, 0) == 3 && A(1, 1) == 4, "TinyMatrix::view_of(data)");

    A.apply(x, y);
    Test(y[0] == 5 && y[1] == 11, "TinyMatrix::apply(vector, vector)");

    TinyMatrix<double, 2, 3> B({ {1,2,3},{4,5,6} });
    TinyMatrix<double, 3, 2> C({ {1,2},{3,4}, {5,6} });
    B.apply(C, A);//A=B*C
    Test(A(0, 0) == 22 && A(0, 1) == 28 && A(1, 0) == 49 && A(1, 1) == 64, "TinyMatrix::apply(matrix, matrix)");

    A = { {2,1},{1,3} };
    inverse(A);
    Test(A(0, 0) == 0.6 && A(0, 1) == -0.2 && A(1, 0) == -0.2 && A(1, 1) == 0.4, "inverse(matrix)");

    auto D = eye<double, 4>();
    inverse(D);
    Test(D(0, 0) == 1 && D(1, 1) == 1 && D(2, 2) == 1 && D(3, 3) == 1, "inverse(matrix)");

    for (auto& dij : D)printf("%lf\n", dij);
}
*/
