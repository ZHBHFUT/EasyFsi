#pragma once
#ifdef __cpp_lib_span
#include <span>

namespace EasyLib {
    inline constexpr size_t dynamic_extent = std::dynamic_extent;
    template<typename T, size_t Extent = dynamic_extent>
    using Span = std::span<T, Extent>;
}

#else
#include <climits> // SIZE_MAX
#include <type_traits>
#include <array>

#include "Assert.hpp"
#include "RandomAccessIterator.hpp"

namespace EasyLib {

    inline constexpr size_t dynamic_extent =  SIZE_MAX;

    template<typename T, size_t Extent = dynamic_extent> class Span;

    template<typename T, size_t Extent>
    class SpanBase
    {
    public:
        constexpr SpanBase(T* base, [[maybe_unused]] size_type count)noexcept
            :_Mydata(base)
        {
            ASSERT(count == Extent);
        }

    protected:
        T* _Mydata{ nullptr };
    };
    template<typename T>
    class SpanBase<T, dynamic_extent>
    {
    public:
        constexpr SpanBase(T* base, size_t count)noexcept :_Mydata(base), _Mysize(count) {}
    protected:
        T*     _Mydata{ nullptr };
        size_t _Mysize{ 0 };
    };

    //! @brief Fixed-length span
    //! @tparam T       element type
    //! @tparam Extent  element number
    template<typename T, size_t Extent>
    class Span : public SpanBase<T, Extent>
    {
    public:
        using base_type        = SpanBase<T, Extent>;
        using element_type     = T;
        using value_type       = std::remove_cv_t<T>;
        using size_type        = size_type;
        using difference_type  = std::ptrdiff_t;
        using pointer          = T*;
        using const_pointer    = const T*;
        using reference        = T&;
        using const_reference  = const T&;
        using iterator         = RandomAccessIterator<T>;
        using reverse_iterator = std::reverse_iterator<iterator>;

        inline static constexpr size_type extent = Extent;

        constexpr Span() noexcept = default;
        constexpr Span(const Span&)noexcept = default;

        constexpr Span(T* base, size_type count)
            :base_type(base, count)
        {
            //if constexpr (Extent != dynamic_extent) {
            //    ASSERT(count == Extent);
            //}
        }
        template<size_type N>
        constexpr Span(T(&arr)[N])noexcept
            : base_type(arr, N)
        {
            static_assert(Extent == dynamic_extent || Extent == N, "length not agree!");
        }
        template<size_type N>
        constexpr Span(std::array<T, N>& arr)noexcept
            : base_type(arr.data(), N)
        {
            static_assert(Extent == dynamic_extent || Extent == N, "length not agree!");
        }
        template<typename T2, size_type N>
        constexpr Span(const std::array<T2, N>& arr)noexcept
            : base_type(arr.data(), N)
        {
            static_assert(std::is_const_v<T> && std::is_same_v<value_type, std::remove_cv_t<T2>>, "invalid type!");
        }
        template<typename T2, size_t N>
        constexpr Span(const Span<T2, N>& s)noexcept
            :base_type(s.data(), N)
        {
            static_assert(Extent == dynamic_extent || Extent == N, "length not agree!");
        }
        template<typename Container>
        constexpr Span(Container&& ctn)noexcept
            :Span(ctn.data(), ctn.size())
        {}

        constexpr pointer   data()const noexcept { return base_type::_Mydata; }
        constexpr size_type size()const noexcept
        {
            if constexpr (Extent == dynamic_extent)
                return base_type::_Mysize;
            else
                return extent;
        }

        constexpr Span& operator = (const Span&)noexcept = default;

        constexpr size_type size_bytes()const noexcept { return sizeof(element_type) * size(); }

        constexpr iterator begin ()const noexcept { return iterator{ base_type::_Mydata }; }
        constexpr iterator end   ()const noexcept { return iterator{ base_type::_Mydata + size() }; }
        constexpr iterator rbegin()const noexcept { return reverse_iterator(end()); }
        constexpr iterator rend  ()const noexcept { return reverse_iterator(begin()); }

        constexpr reference front()const noexcept { return *base_type::_Mydata; }
        constexpr reference back ()const noexcept { return *(base_type::_Mydata + size() - 1); }

        constexpr reference operator[](size_type i)const noexcept { ASSERT(i >= 0 && i < size()); return base_type::_Mydata[i]; }

        constexpr bool empty()const noexcept { return size() == 0; }

        template<size_type Count>
        constexpr Span<T, Count> first()const
        {
            ASSERT(Count <= size());
            return Span<T, Count>(base_type::_Mydata, Count);
        }
        constexpr Span<T, dynamic_extent> first(size_type count)const
        {
            ASSERT(count <= size());
            return Span<T, dynamic_extent>(base_type::_Mydata, count);
        }

        template<size_type Count>
        constexpr Span<T, Count> last()const
        {
            ASSERT(Count <= size());
            return Span<T, Count>(base_type::_Mydata + (size() - Count), Count);
        }
        constexpr Span<T, dynamic_extent> last(size_type count)const
        {
            ASSERT(count <= size());
            return Span<T, dynamic_extent>(base_type::_Mydata + (size() - count), count);
        }

        template<size_type Offset, size_type Count>
        constexpr Span<T, Count> subspan()const
        {
            ASSERT(Offset + Count <= size());
            return Span<T, Count>(base_type::_Mydata + Offset, Count);
        }
        constexpr Span<T, dynamic_extent> subspan(size_type offset, size_type count)const
        {
            ASSERT(offset + count <= size());
            return Span<T, dynamic_extent>(base_type::_Mydata + offset, count);
        }
    };

    template<typename T, size_t N>
    auto as_bytes(Span<T> s)
    {
        if constexpr (N == dynamic_extent)
            return Span<const std::byte>(
                reinterpret_cast<const std::byte*>(s.data()),
                s.size_bytes()
            );
        else
            return Span<const std::byte, sizeof(T)* N>(
                reinterpret_cast<const std::byte*>(s.data()),
                sizeof(T) * N
            );
    }
    template<typename T, size_t N>
    auto as_writable_bytes(Span<T> s)
    {
        static_assert(!std::is_const_v<T>, "type is not writable!");

        if constexpr (N == dynamic_extent)
            return Span<std::byte>(
                reinterpret_cast<std::byte*>(s.data()),
                s.size_bytes()
            );
        else
            return Span<std::byte, sizeof(T)* N>(
                reinterpret_cast<std::byte*>(s.data()),
                sizeof(T) * N
            );
    }
}
#endif

namespace EasyLib {

    template<typename Container>
    Span<typename Container::value_type> make_span(Container& ctn)
    {
        using value_type = typename Container::value_type;
        if constexpr (std::is_const_v<Container>)
            return Span<const value_type>(ctn.data(), ctn.size());
        else
            return Span<value_type>(ctn.data(), ctn.size());
    }

    template<typename T, size_t N>
    Span<T, N> make_span(T(&pointer)[N])
    {
        return Span<T, N>(pointer);
    }

    template<typename T>
    Span<T, dynamic_extent> make_span(T* pointer, size_t count)
    {
        return Span<T, dynamic_extent>{pointer, count};
    }
}
