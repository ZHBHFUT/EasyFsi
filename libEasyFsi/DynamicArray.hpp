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
//! @file       DynamicArray.hpp
//!             The definition of DynamicArray class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <type_traits>
#include <algorithm>
#include <numeric>
#include <array>
#include <stdexcept>
#include <vector>

#include "Assert.hpp"
#include "Index.hpp"
#include "Inline.hpp"
#include "Span.hpp"

namespace EasyLib {

#if __cplusplus < 201703L
    namespace detail {
        _force_inline_ int_l product() { return 1; }
        template<typename ... Args>
        _force_inline_ int_l product(int_l x0, Args ... args) { return x0 * product(args...); }
    }
#endif

    template<typename T, int RANK> class DynamicArray;

    template<typename T, int RANK>
    class DynamicArray
    {
    public:
        using element_type    = T;
        using size_type       = int_l;
        using value_type      = std::remove_cv_t<T>;
        using reference       = T&;
        using pointer         = T*;
        using const_reference = const value_type&;
        using const_pointer   = const value_type*;
        using storage         = std::vector<T>;
        using iterator        = typename storage::iterator;
        using const_iterator  = typename storage::const_iterator;
        using reverse_iterator       = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    
        static_const int Rank = RANK;
    
        DynamicArray() = default;
    
        DynamicArray(const DynamicArray& other) = default;
    
        DynamicArray(DynamicArray&& other)noexcept
            :data_(std::move(other.data_)),
            extent_(std::move(other.extent_)),
            stride_(std::move(other.stride_))
        {}
    
        DynamicArray& operator = (const DynamicArray& other) = default;
    
        DynamicArray& operator = (DynamicArray&& other)noexcept
        {
            data_   = std::move(other.data_);
            extent_ = std::move(other.extent_);
            stride_ = std::move(other.stride_);
            return *this;
        }
    
        ~DynamicArray() = default;
    
        template<typename ... SizeTypes>
        DynamicArray(size_type m, size_type n, SizeTypes ... other_sizes)
            :
#if __cplusplus < 201703L
            data_(detail::product(m, n, other_sizes...)),
#else
            data_(m* (n *...* other_sizes)),
#endif
            extent_{ m,n,other_sizes... }
        {
            static_assert(sizeof...(other_sizes) + 2 == RANK, "invalid argument number");

            // update size info
            update_stride_();
            ASSERT(extent_[0] * stride_[0] == data_.size());
        }
    
        int           rank  ()const noexcept { return RANK; }
        size_type     extent(int rank_id)const  { return extent_[rank_id]; }
        size_type     numel ()const noexcept { return data_.size(); }
        pointer       data  ()      noexcept { return data_.data(); }
        const_pointer data  ()const noexcept { return data_.data(); }
        bool          empty ()const noexcept { return data_.empty(); }
    
        void clear()noexcept
        {
            data_.clear();
            extent_.fill(0);
            stride_.fill(0);
        }
    
        void reserve(size_type n)
        {
            data_.reserve(n);
        }

        template<typename ... SizeTypes>
        void reserve(size_type m, size_type n, SizeTypes ... other_sizes)
        {
            static_assert(sizeof...(other_sizes) + 2 == RANK, "invalid argument number");
            data_.reserve(
#if __cplusplus < 201703L
                detail::product(m, n, other_sizes...)
#else
                m * (n * ... * other_sizes)
#endif
            );
        }

        template<typename ... SizeTypes>
        void resize(size_type m, size_type n, SizeTypes ... other_sizes)
        {
            static_assert(sizeof...(other_sizes) + 2 == RANK, "invalid argument number");

            data_.resize(
#if __cplusplus < 201703L
                detail::product(m, n, other_sizes...)
#else
                m * (n * ... * other_sizes)
#endif
            );
    
            // update size info
            extent_ = { m,n,other_sizes... };
            update_stride_();
            ASSERT(extent_[0] * stride_[0] == data_.size());
        }
    
        void resize(const std::array<size_type, RANK>& size)
        {
            // resizing
            data_.resize(std::reduce(size.begin(), size.end(), size_type{ 1 }, std::multiplies<size_type>()));
            
            // update size info
            extent_ = size;
            update_stride_();
            ASSERT(extent_[0] * stride_[0] == data_.size());
        }
    
        template<typename ... Args>
        void reshape(size_type m, size_type n, Args ... other_sizes)
        {
            static_assert(sizeof...(other_sizes) + 2 == RANK, "invalid argument number");
            auto nelem = 
#if __cplusplus < 201703L
                detail::product(m, n, other_sizes...)
#else
                m * (n *...* other_sizes)
#endif
                ;
            if (nelem != data_.size())throw std::invalid_argument("DynamicArray::reshape(), size not agree!");
    
            extent_ = { m,n,other_sizes... };
            update_stride_();
            ASSERT(extent_[0] * stride_[0] == data_.size());
        }
    
        void reshape(const std::array<size_type, RANK>& size)
        {
            size_type nelem = std::reduce(size.begin(), size.end(), size_type{ 1 }, std::multiplies<size_type>());
            if (nelem != data_.size())throw std::invalid_argument("DynamicArray::reshape(), size not agree!");
    
            extent_ = size;
            update_stride_();
            ASSERT(extent_[0] * stride_[0] == data_.size());
        }
    
        void fill(const value_type& value)noexcept
        {
            std::fill(data_.begin(), data_.end(), value);
        }
    
        void swap(DynamicArray& other)noexcept
        {
            std::swap(data_, other.data_);
            std::swap(extent_, other.extent_);
            std::swap(stride_, other.stride_);
        }
    
        void swap_elements(DynamicArray& other)
        {
            if (data_.size() != other.data_.size())
                throw std::logic_error("DynamicArray::swap_elements, size not agree!");
            std::swap_ranges(data_.begin(), data_.end(), other.data_.begin());
        }

        void copy_elements(const DynamicArray& other)
        {
            if (numel() != other.numel())throw std::logic_error("DynamicArray::copy_elements, size not agree!");
            std::copy(other.begin(), other.end(), begin());
        }
    
        reference       operator()(size_type i)      noexcept { return data_[i]; }
        const_reference operator()(size_type i)const noexcept { return data_[i]; }
    
        reference       operator()(const std::array<size_type, RANK>& idx)      noexcept { return data_[map_1d(idx)]; }
        const_reference operator()(const std::array<size_type, RANK>& idx)const noexcept { return data_[map_1d(idx)]; }
        reference       operator[](const std::array<size_type, RANK>& idx)      noexcept { return data_[map_1d(idx)]; }
        const_reference operator[](const std::array<size_type, RANK>& idx)const noexcept { return data_[map_1d(idx)]; }

        template<typename ... Args>
        reference       operator()(size_type i, size_type j, Args ... args)noexcept
        {
            return data_[map_1d(i, j, args...)];
        }
        template<typename ... Args>
        const_reference operator()(size_type i, size_type j, Args ... args)const noexcept
        {
            return data_[map_1d(i, j, args...)];
        }
    
        template<typename ... Args>
        size_type map_1d(size_type i0, size_type i1, Args ... args)const noexcept
        {
            static_assert(sizeof...(args) + 2 == RANK, "invalid arguments number!");
            return map_1d_imp_(i0, i1, args...);
        }
        size_type map_1d(const std::array<size_type, RANK>& idx)const noexcept
        {
            return std::inner_product(idx.begin(), idx.end(), stride_.begin(), size_type{ 0 });
        }
    
        iterator       begin()noexcept { return data_.begin(); }
        iterator       end  ()noexcept { return data_.end(); }
        const_iterator begin()const noexcept { return data_.begin(); }
        const_iterator end  ()const noexcept { return  data_.end(); }
        reverse_iterator       rbegin()noexcept { return reverse_iterator(end()); }
        reverse_iterator       rend  ()noexcept { return reverse_iterator(begin()); }
        const_reverse_iterator rbegin()const noexcept { return const_reverse_iterator(end()); }
        const_reverse_iterator rend  ()const noexcept { return const_reverse_iterator(begin()); }
    
        template<typename, int>
        friend class DynamicArray;
        friend class Communicator;

    private:
        template<typename ... Args>
        _force_inline_ size_type map_1d_imp_(size_type i0, Args ... args)const noexcept
        {
            static_assert(sizeof...(args) > 0, "");
            constexpr size_type idim = RANK - 1 - sizeof...(args);
            ASSERT(i0 >= 0 && i0 < std::get<idim>(extent_));
            return i0 * std::get<idim>(stride_) + map_1d_imp_(args...);
        }
        _force_inline_ size_type map_1d_imp_(size_type i0)const noexcept
        {
            static constexpr size_type idim = RANK - 1;
            ASSERT(i0 >= 0 && i0 < std::get<idim>(extent_));
            return i0 * stride_.back();
        }
        void update_stride_()noexcept
        {
            stride_.back() = 1;
            for (int idim = 1; idim < RANK; ++idim) {
                stride_[RANK - idim - 1] = extent_[RANK - idim] * stride_[RANK - idim];
            }
        }

    private:
        std::vector<T>              data_;
        std::array<size_type, RANK> extent_  { 0 };
        std::array<size_type, RANK> stride_{ 0 };
    };

    template<typename T>
    class DynamicArray<T, 1> : public std::vector<T>
    {
    public:
        using storage         = std::vector<T>;
        using element_type    = T;
        using value_type      = typename storage::value_type;
        using size_type       = typename storage::size_type;
        using pointer         = T*;
        using reference       = T&;
        using const_pointer   = const value_type*;
        using const_reference = const value_type&;

        static_const int Rank = 1;
        
        using storage::storage;
        using storage::operator=;
        using storage::operator[];
        using storage::size;
        using storage::empty;
        using storage::begin;
        using storage::end;
        using storage::resize;

        constexpr int           rank  ()const noexcept { return 1; }
        size_type     extent(int idim = 0)const { (void)idim; return size(); }
        size_type     numel ()const noexcept { return size(); }

        void fill(const value_type& value)
        {
            std::fill(begin(), end(), value);
        }

        void swap_elements(DynamicArray& other)
        {
            if (size() != other.size())throw std::logic_error("DynamicArray::swap_elements, size not agree!");
            std::swap_ranges(begin(), end(), other.begin());
        }

        void copy_elements(const DynamicArray& other)
        {
            if (size() != other.size())throw std::logic_error("DynamicArray::copy_elements, size not agree!");
            std::copy(other.begin(), other.end(), begin());
        }

        reference       operator()(size_type i)      noexcept { return operator[](i); }
        const_reference operator()(size_type i)const noexcept { return operator[](i); }

        template<typename, int>
        friend class DynamicArray;
        friend class Communicator;
    };
    
    template<typename T>
    class DynamicArray<T, 2>
    {
    public:
        using storage         = std::vector<T>;
        using element_type    = T;
        using value_type      = typename storage::value_type;
        using size_type       = typename storage::size_type;
        using pointer         = T*;
        using reference       = T&;
        using const_pointer   = const value_type*;
        using const_reference = const value_type&;
        using iterator        = typename storage::iterator;
        using const_iterator  = typename storage::const_iterator;
        using reverse_iterator = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;

        static_const int Rank = 2;

        DynamicArray() = default;
        DynamicArray(const DynamicArray&) = default;
        DynamicArray(DynamicArray&& other)noexcept
            :data_(std::move(other.data_)),
            extent_(std::move(other.extent_))
        {}
        DynamicArray(size_type m, size_type n, const_reference value = value_type{})
            :data_(m * n, value),
            extent_{m, n}
        {}

        DynamicArray& operator = (const DynamicArray&) = default;
        DynamicArray& operator = (DynamicArray&& other)noexcept
        {
            data_           = std::move(other.data_);
            extent_         = std::move(other.extent_);
            return *this;
        }

        ~DynamicArray() = default;

        int           rank()const noexcept { return 2; }
        size_type     extent(int idim)const { ASSERT(idim >= 0 && idim < Rank); return extent_[idim]; }
        size_type     numel()const noexcept { return extent_[0] * extent_[1]; }
        pointer       data()      noexcept { return data_.data(); }
        const_pointer data()const noexcept { return data_.data(); }
        bool          empty()const noexcept { return extent_[0] == 0 || extent_[1] == 0; }

        void clear()noexcept
        {
            data_.clear();
            extent_[0] = extent_[1] = 0;
        }

        void reserve(size_type n)
        {
            data_.reserve(n);
        }

        void resize(size_type m, size_type n, const_reference value = value_type{})
        {
            data_.resize(m * n, value);
            extent_[0] = m;
            extent_[1] = n;
        }

        void fill(const value_type& value)
        {
            std::fill(data_.begin(), data_.end(), value);
        }

        void swap(DynamicArray& other)noexcept
        {
            data_.swap(other.data_);
            std::swap(extent_, other.extent_);
        }

        void swap_elements(DynamicArray& other)
        {
            if (numel() != other.numel())throw std::logic_error("DynamicArray::swap_elements, size not agree!");
            std::swap_ranges(data_.begin(), data_.end(), other.data_.begin());
        }

        void copy_elements(const DynamicArray& other)
        {
            if (numel() != other.numel())throw std::logic_error("DynamicArray::copy_elements, size not agree!");
            std::copy(other.begin(), other.end(), begin());
        }

        reference       operator()(size_type i)      noexcept { ASSERT(i >= 0 && i < extent_[0] * extent_[1]); return data_[i]; }
        const_reference operator()(size_type i)const noexcept { ASSERT(i >= 0 && i < extent_[0] * extent_[1]); return data_[i]; }
        auto operator[](size_type i)      noexcept { ASSERT(i >= 0 && i < extent_[0] * extent_[1]); return make_span(data() + i * extent_[1], extent_[1]); }
        auto operator[](size_type i)const noexcept { ASSERT(i >= 0 && i < extent_[0] * extent_[1]); return make_span(data() + i * extent_[1], extent_[1]); }
        reference       operator()(size_type i, size_type j)      noexcept { ASSERT(i >= 0 && i < extent_[0] && j >= 0 && j < extent_[1]); return data_[i * extent_[1] + j]; }
        const_reference operator()(size_type i, size_type j)const noexcept { ASSERT(i >= 0 && i < extent_[0] && j >= 0 && j < extent_[1]); return data_[i * extent_[1] + j]; }
        //reference       operator[](size_type i, size_type j)      noexcept { ASSERT(i >= 0 && i < extent_[0] && j >= 0 && j < extent_[1]); return data_[i * extent_[1] + j]; }
        //const_reference operator[](size_type i, size_type j)const noexcept { ASSERT(i >= 0 && i < extent_[0] && j >= 0 && j < extent_[1]); return data_[i * extent_[1] + j]; }

        reference       at(size_type i)
        {
            return data_.at(i);
        }
        const_reference at(size_type i)const
        {
            return data_.at(i);
        }
        reference       at(size_type i, size_type j)
        {
            if (i < 0 || i >= extent_[0] ||
                j < 0 || j >= extent_[1])
                throw std::out_of_range("DynamicArray::at(), index out of range!");
            return data_[i * extent_[1] + j];
        }
        const_reference at(size_type i, size_type j)const
        {
            if (i < 0 || i >= extent_[0] ||
                j < 0 || j >= extent_[1])
                throw std::out_of_range("DynamicArray::at(), index out of range!");
            return data_[i * extent_[1] + j];
        }

        size_type map_1d(size_type i, size_type j)const noexcept
        {
            ASSERT(i >= 0 && i < extent_[0] && j >= 0 && j < extent_[1]);
            return i * extent_[1] + j;
        }
        size_type map_1d(const std::array<size_type, 2>& idx)const noexcept
        {
            ASSERT(idx[0] >= 0 && idx[0] < extent_[0] && idx[1] >= 0 && idx[1] < extent_[1]);
            return idx[0] * extent_[1] + idx[1];
        }

        iterator       begin()noexcept { return data_.begin(); }
        iterator       end  ()noexcept { return data_.end(); }
        const_iterator begin()const noexcept { return data_.begin(); }
        const_iterator end  ()const noexcept { return data_.end(); }
        reverse_iterator       rbegin()noexcept { return reverse_iterator(end()); }
        reverse_iterator       rend  ()noexcept { return reverse_iterator(begin()); }
        const_reverse_iterator rbegin()const noexcept { return const_reverse_iterator(end()); }
        const_reverse_iterator rend  ()const noexcept { return const_reverse_iterator(begin()); }

        template<typename, int>
        friend class DynamicArray;
        friend class Communicator;

    private:
        storage                  data_;
        std::array<size_type, 2> extent_{ 0 };
    };
}
