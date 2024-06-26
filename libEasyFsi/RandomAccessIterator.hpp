#pragma once
#include <type_traits> // remove_const_t
#include <iterator>    // random_access_iterator_tag

namespace EasyLib {

    //! @brief A random access iterator class.
    //! @tparam T The element type.
    template<typename T> struct RandomAccessIterator //: std::iterator< std::random_access_iterator_tag, T>
    {
#ifdef _MSC_VER 
        using iterator_category = std::random_access_iterator_tag;
#endif
        using value_type        = std::remove_cv_t<T>;
        using difference_type   = ptrdiff_t;
        using pointer           = T*;
        using reference         = T&;

        RandomAccessIterator() = default;
        RandomAccessIterator(const RandomAccessIterator&) = default;
        RandomAccessIterator& operator=(const RandomAccessIterator&) = default;

        const RandomAccessIterator(T* data, size_t offset = 0)noexcept :ptr_(data ? data + offset : nullptr) {}

        const reference operator *()const noexcept { return *ptr_; }
        const pointer   operator->()const noexcept { return  ptr_; }

        const reference operator [](size_t i)const noexcept { return ptr_[i]; }

        const bool operator ==(const RandomAccessIterator& it)const noexcept { return ptr_ == it.ptr_; }
        const bool operator !=(const RandomAccessIterator& it)const noexcept { return ptr_ != it.ptr_; }
        const bool operator < (const RandomAccessIterator& it)const noexcept { return ptr_ <  it.ptr_; }
        const bool operator <=(const RandomAccessIterator& it)const noexcept { return ptr_ <= it.ptr_; }
        const bool operator > (const RandomAccessIterator& it)const noexcept { return ptr_ >  it.ptr_; }
        const bool operator >=(const RandomAccessIterator& it)const noexcept { return ptr_ >= it.ptr_; }

        const RandomAccessIterator& operator++()noexcept { ++ptr_; return *this; }
        const RandomAccessIterator& operator--()noexcept { --ptr_; return *this; }
        const RandomAccessIterator  operator++(int)noexcept { return RandomAccessIterator{ ptr_ + 1 }; }
        const RandomAccessIterator  operator--(int)noexcept { return RandomAccessIterator{ ptr_ - 1 }; }
        const RandomAccessIterator& operator+=(difference_type _Off)noexcept { ptr_ += _Off; return *this; }
        const RandomAccessIterator& operator-=(difference_type _Off)noexcept { ptr_ -= _Off; return *this; }
        const RandomAccessIterator operator + (difference_type _Off)const noexcept { return RandomAccessIterator{ ptr_ + _Off }; }
        const RandomAccessIterator operator - (difference_type _Off)const noexcept { return RandomAccessIterator{ ptr_ - _Off }; }

        const difference_type operator -(const RandomAccessIterator& it)const noexcept { return static_cast<difference_type>(ptr_ - it.ptr_); }

    private:
        T* ptr_{ nullptr };
    };
} // EasyLib
