#include "LinAlgs.h"
#include "DynamicMatrix.hpp"

namespace EasyLib {

    DynamicMatrix::DynamicMatrix(size_type m, size_type n, value_type init_val/* = value_type{ 0 }*/)
        :size_{ m >= 0 ? m : 0,n >= 0 ? n : 0 },
        data_(static_cast<size_t>(size_[0] * size_[1]), init_val)
    {}

    DynamicMatrix::DynamicMatrix(DynamicMatrix&& m)noexcept
        :size_(m.size_),
        data_(std::move(m.data_))
    {}

    DynamicMatrix& DynamicMatrix::operator=(DynamicMatrix&& m)noexcept
    {
        if (this != std::addressof(m)) {
            size_ = m.size_;
            data_ = std::move(m.data_);
        }
        return *this;
    }

    void DynamicMatrix::clear()noexcept
    {
        size_[0] = size_[1] = 0;
        data_.clear();
    }

    void DynamicMatrix::resize(size_type m, size_type n, value_type init_val/* = value_type{ 0 }*/)
    {
        size_[0] = m >= 0 ? m : 0;
        size_[1] = n >= 0 ? n : 0;
        data_.resize(size_[0] * size_[1], init_val);
    }

    void DynamicMatrix::reshape(size_type m, size_type n)
    {
        m = m >= 0 ? m : 0;
        n = n >= 0 ? n : 0;

        if (m * n != size_[0] * size_[1])
            throw std::invalid_argument("size not agree size");
        size_[0] = m; size_[1] = n;
    }

    void DynamicMatrix::fill(value_type value)
    {
        std::fill(data_.begin(), data_.end(), value);
    }

    void DynamicMatrix::copy_from(const DynamicMatrix& mat)
    {
        if (numel() != mat.numel())
            resize(mat.size_[0], mat.size_[1]);
        else {
            size_[0] = mat.size_[0];
            size_[1] = mat.size_[1];
        }
        std::copy(mat.begin(), mat.end(), begin());
    }

    void DynamicMatrix::identity()
    {
        fill(0);
        for (size_type i = 0; i < std::min(size_[0], size_[1]); ++i)
            data_[i * size_[1] + i] = 1;
    }
    bool DynamicMatrix::inverse()
    {
        if (size_[0] != size_[1])throw std::runtime_error("matrix is not square!");

        int singular = 0;
        std::vector<int> buffer(2 * size_[0]);
        mat_inverse((int)size_[0], data_.data(), buffer.data(), &singular);
        return singular == 0;
    }
    bool DynamicMatrix::inverse(std::vector<int>& buffer)
    {
        if (size_[0] != size_[1])throw std::runtime_error("matrix is not square!");

        int singular = 0;
        buffer.resize(2 * size_[0]);
        mat_inverse((int)size_[0], data_.data(), buffer.data(), &singular);
        return singular == 0;
    }
    //! @brief compute y = A.x
    void DynamicMatrix::apply(const DynamicVector& x, DynamicVector& y)const
    {
        // check size
        if (x.size() != size_[1] ||
            y.size() != size_[0]) {
            ASSERT(false);
            throw std::out_of_range("dimension not agree.");
        }
        mat_apply_vec((int)size_[0], (int)size_[1], data_.data(), x.data(), y.data());
    }
    void DynamicMatrix::apply(const DynamicMatrix& x, DynamicMatrix& y)const
    {
        // check size
        if (x.size(0) != size_[1] ||
            x.size(1) != y.size(1) ||
            y.size(0) != size_[0]) {
            ASSERT(false);
            throw std::out_of_range("dimension not agree.");
        }
        mat_apply_mat((int)size_[0], (int)size_[1], (int)x.size(1), data_.data(), x.data(), y.data());
    }

    //! @brief compute y = A . x + y
    void DynamicMatrix::apply_add(const DynamicVector& x, DynamicVector& y)const
    {
        // check size
        if (x.size() != size_[1] ||
            y.size() != size_[0]) {
            ASSERT(false);
            throw std::out_of_range("dimension not agree.");
        }
        mat_apply_add_vec((int)size_[0], (int)size_[1], data_.data(), x.data(), y.data());
    }
    void DynamicMatrix::apply_add(const DynamicMatrix& x, DynamicMatrix& y)const
    {
        // check size
        if (x.size(0) != size_[1] ||
            x.size(1) != y.size(1) ||
            y.size(0) != size_[0]) {
            ASSERT(false);
            throw std::out_of_range("dimension not agree.");
        }
        mat_apply_add_mat((int)size_[0], (int)size_[1], (int)x.size(1), data_.data(), x.data(), y.data());
    }

}
