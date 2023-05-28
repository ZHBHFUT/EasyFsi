#pragma once
#include <type_traits>
#include <iterator>

namespace EasyLib {
    template<typename T> struct RandomAccessIterator //: std::iterator< std::random_access_iterator_tag, T>
    {
        using iterator_category = std::random_access_iterator_tag;
        using value_type = T;
        using difference_type = ptrdiff_t;
        using pointer = T*;
        using reference = T&;
        using const_reference = const std::remove_const_t<T>&;

        RandomAccessIterator(T* data)noexcept :ptr_(data) {}
        RandomAccessIterator(const RandomAccessIterator&) = default;
        RandomAccessIterator& operator=(const RandomAccessIterator&) = default;

        inline bool operator ==(const RandomAccessIterator& it)const noexcept { return ptr_ == it.ptr_; }
        inline bool operator !=(const RandomAccessIterator& it)const noexcept { return ptr_ != it.ptr_; }
        inline bool operator < (const RandomAccessIterator& it)const noexcept { return ptr_ < it.ptr_; }
        inline bool operator <=(const RandomAccessIterator& it)const noexcept { return ptr_ <= it.ptr_; }
        inline bool operator > (const RandomAccessIterator& it)const noexcept { return ptr_ > it.ptr_; }
        inline bool operator >=(const RandomAccessIterator& it)const noexcept { return ptr_ >= it.ptr_; }

        inline RandomAccessIterator& operator++()noexcept { ++ptr_; return *this; }
        inline RandomAccessIterator& operator--()noexcept { --ptr_; return *this; }
        inline RandomAccessIterator  operator++(int)noexcept { return RandomAccessIterator{ ptr_ + 1 }; }
        inline RandomAccessIterator  operator--(int)noexcept { return RandomAccessIterator{ ptr_ - 1 }; }
        inline RandomAccessIterator& operator+=(difference_type _Off)noexcept { ptr_ += _Off; return *this; }
        inline RandomAccessIterator& operator-=(difference_type _Off)noexcept { ptr_ -= _Off; return *this; }
        inline RandomAccessIterator operator + (difference_type _Off)noexcept { return RandomAccessIterator{ ptr_ + _Off }; }
        inline RandomAccessIterator operator - (difference_type _Off)noexcept { return RandomAccessIterator{ ptr_ - _Off }; }

        inline ptrdiff_t operator -(const RandomAccessIterator& it)const noexcept { return static_cast<ptrdiff_t>(ptr_ - it.ptr_); }

        inline reference operator *()      noexcept { return *ptr_; }
        inline pointer   operator->()const noexcept { return  ptr_; }

        inline       reference operator [](size_t i) { return ptr_[i]; }
        inline const_reference operator [](size_t i)const { return ptr_[i]; }

    private:
        T* ptr_{ nullptr };
    };
}
