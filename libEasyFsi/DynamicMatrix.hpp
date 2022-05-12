#pragma once
#include <algorithm>
#include <stdexcept>
#include <array>
#include <vector>
#include <span>

#include "Assert.hpp"
#include "DynamicVector.hpp"

namespace EasyLib {

    class Communicator;

    class DynamicMatrix
    {
    public:
        using value_type             = double;
        using size_type              = size_t;
        using reference              = double&;
        using const_reference        = const double&;
        using pointer                = double*;
        using const_pointer          = const double*;

        using iterator               = std::vector<double>::iterator;
        using const_iterator         = std::vector<double>::const_iterator;
        using reverse_iterator       = std::vector<double>::reverse_iterator;
        using const_reverse_iterator = std::vector<double>::const_reverse_iterator;

        inline static constexpr int dimension = 2;

        DynamicMatrix() = default;
        DynamicMatrix(const DynamicMatrix&) = default;
        DynamicMatrix& operator=(const DynamicMatrix&) = default;

        DynamicMatrix(DynamicMatrix&& m)noexcept;

        DynamicMatrix& operator=(DynamicMatrix&& m)noexcept;

        DynamicMatrix(size_type m, size_type n, value_type init_val = value_type{ 0 });

        void clear()noexcept;

        void resize(size_type m, size_type n, value_type init_val = value_type{ 0 });

        void reshape(size_type m, size_type n);

        void fill(value_type value);

        void copy_from(const DynamicMatrix& mat);

        inline constexpr int  ndim()const noexcept { return 2; }

        inline constexpr bool is_square()const noexcept { return size_[0] == size_[1]; }

        inline constexpr bool empty()const noexcept { return data_.empty(); }

        inline size_type nrow()const noexcept { return size_[0]; }
        inline size_type ncol()const noexcept { return size_[1]; }

        inline size_type size(int idim)const noexcept { ASSERT(idim >= 0 && idim < 2); return size_[idim]; }
        inline size_type numel()const noexcept { return size_[0] * size_[1]; }

        inline       pointer data()      noexcept { return data_.data(); }
        inline const_pointer data()const noexcept { return data_.data(); }

        inline auto operator[](size_type i)
        {
            ASSERT(i >= 0 && i < size_[0]);
            return std::span<value_type>{&data_.at(i* size_[1]), size_[1]};
        }
        inline auto operator[](size_type i)const
        {
            ASSERT(i >= 0 && i < size_[0]);
            return std::span<const value_type>{&data_.at(i* size_[1]), size_[1]};
        }

        inline reference       operator()(size_type i)                   noexcept { ASSERT(i >= 0 && i < data_.size()); return data_[i]; }
        inline const_reference operator()(size_type i)const              noexcept { ASSERT(i >= 0 && i < data_.size()); return data_[i]; }
        inline reference       operator()(size_type i, size_type j)      noexcept { ASSERT(i >= 0 && i < size_[0] && j >= 0 && j < size_[1]); return data_[i * size_[1] + j]; }
        inline const_reference operator()(size_type i, size_type j)const noexcept { ASSERT(i >= 0 && i < size_[0] && j >= 0 && j < size_[1]); return data_[i * size_[1] + j]; }

        inline reference       at(size_type i) { return data_.at(i); }
        inline const_reference at(size_type i)const { return data_.at(i); }
        inline reference       at(size_type i, size_type j) { ASSERT(i >= 0 && i < size_[0] && j >= 0 && j < size_[1]); return data_.at(i * size_[1] + j); }
        inline const_reference at(size_type i, size_type j)const { ASSERT(i >= 0 && i < size_[0] && j >= 0 && j < size_[1]); return data_.at(i * size_[1] + j); }

        inline auto begin ()noexcept { return data_.begin(); }
        inline auto end   ()noexcept { return data_.end(); }
        inline auto begin ()const noexcept { return data_.begin(); }
        inline auto end   ()const noexcept { return data_.end(); }
        inline auto cbegin()const noexcept { return data_.cbegin(); }
        inline auto cend  ()const noexcept { return data_.cend(); }

        inline auto rbegin()noexcept { return data_.rbegin(); }
        inline auto rend  ()noexcept { return data_.rend(); }
        inline auto rbegin()const noexcept { return data_.rbegin(); }
        inline auto rend  ()const noexcept { return data_.rend(); }

        //! 设置为单位矩阵
        //! @note 如果矩阵不是方阵，则填充对角线元素为1
        void identity();

        //! @brief 对矩阵求逆
        //! @return true=成功，false=失败（矩阵奇异）
        bool inverse();
        bool inverse(std::vector<int>& buffer);

        //! @brief compute y = A.x
        void apply(const DynamicVector& x, DynamicVector& y)const;
        void apply(const DynamicMatrix& x, DynamicMatrix& y)const;

        //! @brief compute y = A . x + y
        void apply_add(const DynamicVector& x, DynamicVector& y)const;
        void apply_add(const DynamicMatrix& x, DynamicMatrix& y)const;

        friend class Communicator;

    private:
        std::array<size_type,2> size_{ 0,0 };
        std::vector<value_type> data_;
    };

}
