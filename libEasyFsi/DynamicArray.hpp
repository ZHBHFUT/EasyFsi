#pragma once
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <span>
#include <array>
#include <stdexcept>
#include <vector>

#include "Assert.hpp"
#include "Index.hpp"
#include "Inline.hpp"

namespace EasyLib {

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
    
        inline static constexpr int Rank = RANK;
    
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
            data_(m * (n *...* other_sizes)),
            extent_{ m,n,other_sizes... }
        {
            static_assert(sizeof...(other_sizes) + 2 == RANK, "invalid argument number");
    
            // update size info
            stride_[RANK - 1] = 1;
            for (size_type i = 2; i <= RANK; ++i)
                stride_[RANK - i] = extent_[RANK - i] * stride_[RANK - i + 1];
            }
    
        inline constexpr int           rank  ()const noexcept { return RANK; }
        inline constexpr size_type     extent(int rank_id)const  { return extent_[rank_id]; }
        inline constexpr size_type     numel ()const noexcept { return data_.size(); }
        inline constexpr pointer       data  ()      noexcept { return data_.data(); }
        inline constexpr const_pointer data  ()const noexcept { return data_.data(); }
        inline constexpr bool          empty ()const noexcept { return data_.empty(); }
    
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
            data_.reserve(m * (n * ... * other_sizes));
        }

        template<typename ... SizeTypes>
        void resize(size_type m, size_type n, SizeTypes ... other_sizes)
        {
            static_assert(sizeof...(other_sizes) + 2 == RANK, "invalid argument number");

            data_.resize(m * (n * ... * other_sizes));
    
            // update size info
            extent_ = { m,n,other_sizes... };
            stride_[RANK - 1] = 1;
            for (int i = 2; i <= RANK; ++i)
                stride_[RANK - i] = extent_[RANK - i] * stride_[RANK - i + 1];
        }
    
        void resize(const std::array<size_type, RANK>& size)
        {
            // resizing
            data_.resize(std::reduce(size.begin(), size.end(), size_type{ 1 }, std::multiplies<size_type>()));
            
            // update size info
            extent_ = size;
            stride_[RANK - 1] = 1;
            for (int i = 2; i <= RANK; ++i)
                stride_[RANK - i] = extent_[RANK - i] * stride_[RANK - i + 1];
        }
    
        template<typename ... Args>
        void reshape(size_type m, size_type n, Args ... other_sizes)
        {
            static_assert(sizeof...(other_sizes) + 2 == RANK, "invalid argument number");
            auto nelem = m * (n *...* other_sizes);
            if (nelem != data_.size())throw std::invalid_argument("DynamicArray::reshape(), size not agree!");
    
            extent_ = { m,n,other_sizes... };
            stride_[RANK - 1] = 1;
            for (int i = 2; i <= RANK; ++i)
                stride_[RANK - i] = extent_[RANK - i] * stride_[RANK - i + 1];
        }
    
        void reshape(const std::array<size_type, RANK>& size)
        {
            size_type nelem = std::reduce(size.begin(), size.end(), size_type{ 1 }, std::multiplies<size_type>());
            if (nelem != data_.size())throw std::invalid_argument("DynamicArray::reshape(), size not agree!");
    
            extent_ = size;
            stride_[RANK - 1] = 1;
            for (int i = 2; i <= RANK; ++i)
                stride_[RANK - i] = extent_[RANK - i] * stride_[RANK - i + 1];
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
            constexpr size_type idim = RANK - 1 - sizeof...(args);
            ASSERT(i0 >= 0 && i0 < std::get<idim>(extent_));
            if constexpr (sizeof...(Args) > 0)
                return i0 * std::get<idim>(stride_) * map_1d_imp_(args...);
            else
                return i0 * stride_.back();
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

        inline static constexpr int Rank = 1;
        
        using storage::storage;
        using storage::operator=;
        using storage::operator[];
        using storage::size;
        using storage::empty;
        using storage::begin;
        using storage::end;
        using storage::resize;

        inline constexpr int           rank  ()const noexcept { return 1; }
        inline constexpr size_type     extent(int idim = 0)const { (void)idim; return size(); }
        inline constexpr size_type     numel ()const noexcept { return size(); }

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

        inline static constexpr int Rank = 2;

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

        inline constexpr int           rank  ()const noexcept { return 2; }
        inline constexpr size_type     extent(int idim)const  { ASSERT(idim >= 0 && idim < Rank); return extent_[idim]; }
        inline constexpr size_type     numel ()const noexcept { return extent_[0] * extent_[1]; }
        inline constexpr pointer       data  ()      noexcept { return data_.data(); }
        inline constexpr const_pointer data  ()const noexcept { return data_.data(); }
        inline constexpr bool          empty ()const noexcept { return extent_[0] == 0 || extent_[1] == 0; }

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
        auto operator[](size_type i)      noexcept { ASSERT(i >= 0 && i < extent_[0] * extent_[1]); return std::span<    element_type>(data() + i * extent_[1], extent_[1]); }
        auto operator[](size_type i)const noexcept { ASSERT(i >= 0 && i < extent_[0] * extent_[1]); return std::span<const value_type>(data() + i * extent_[1], extent_[1]); }
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
