#pragma once
#ifdef __cpp_lib_span
#include <span>

namespace EasyLib {
    inline const size_t dynamic_extent = std::dynamic_extent;
    template<typename T, size_t Extent = dynamic_extent>
    using Span = std::span<T, Extent>;
}

#else
#include <climits> // SIZE_MAX
#include <type_traits>
#include <array>

#include "Assert.hpp"
#include "Inline.hpp"
#include "RandomAccessIterator.hpp"

namespace EasyLib {

    inline_const std::size_t dynamic_extent =  SIZE_MAX;

    template<typename T, typename = void>
    struct has_data_function : std::false_type {};
    template<typename T>
    struct has_data_function<T, std::void_t<decltype(std::declval<T>().data())>> : std::true_type {};
    template<typename T, typename = void>
    struct has_size_function : std::false_type {};
    template<typename T>
    struct has_size_function<T, std::void_t<decltype(std::declval<T>().size())>> : std::true_type {};
    template<typename T>
    using has_data_size_functions = std::conjunction<has_data_function<T>, has_size_function<T>>;

    template<typename T, std::size_t Extent = dynamic_extent> class Span;

    template<typename T, std::size_t Extent>
    class SpanBase
    {
    public:
        using size_type = std::size_t;
        SpanBase(T* base, size_type count)noexcept
            :_Mydata(base)
        {
            ASSERT(count == Extent);
        }
        const size_type size()const noexcept { return Extent; }
    protected:
        T* _Mydata{ nullptr };
    };
    template<typename T>
    class SpanBase<T, dynamic_extent>
    {
    public:
        using size_type = std::size_t;
        SpanBase()noexcept = default;
        SpanBase(T* base, size_type count)noexcept :_Mydata(base), _Mysize(count) {}
        size_type size()const noexcept { return _Mysize; }
    protected:
        T*        _Mydata{ nullptr };
        size_type _Mysize{ 0 };
    };

    //! @brief Fixed-length span
    //! @tparam T       element type
    //! @tparam Extent  element number
    template<typename T, std::size_t Extent>
    class Span : public SpanBase<T, Extent>
    {
    public:
        using base_type        = SpanBase<T, Extent>;
        using element_type     = T;
        using value_type       = std::remove_cv_t<T>;
        using size_type        = base_type::size_type;
        using difference_type  = std::ptrdiff_t;
        using pointer          = T*;
        using const_pointer    = const T*;
        using reference        = T&;
        using const_reference  = const T&;
        using iterator         = RandomAccessIterator<T>;
        using reverse_iterator = std::reverse_iterator<iterator>;

        using base_type::size;

        static_const size_type extent = Extent;

        Span()noexcept = default;
        Span(const Span&)noexcept = default;

        Span(T* base, size_type count)noexcept
            :base_type(base, count)
        {}
        template<size_type N>
        Span(T(&arr)[N])noexcept
            : base_type(arr, N)
        {
            static_assert(Extent == dynamic_extent || Extent == N, "length not agree!");
        }
        template<size_type N>
        Span(std::array<T, N>& arr)noexcept
            : base_type(arr.data(), N)
        {
            static_assert(Extent == dynamic_extent || Extent == N, "length not agree!");
        }
        template<typename T2, size_type N>
        Span(const std::array<T2, N>& arr)noexcept
            : base_type(arr.data(), N)
        {
            static_assert(std::is_const_v<T> && std::is_same_v<value_type, std::remove_cv_t<T2>>, "invalid type!");
        }
        template<typename T2, std::size_t N>
        Span(const Span<T2, N>& s)noexcept
            :base_type(s.data(), N)
        {
            static_assert(Extent == dynamic_extent || Extent == N, "length not agree!");
        }
        template<typename Container, typename = std::enable_if_t<has_data_size_function<Container>::value>>
        Span(Container&& ctn)noexcept
            :Span(ctn.data(), ctn.size())
        {}

        const pointer   data()const noexcept { return base_type::_Mydata; }
        
        Span& operator = (const Span&)noexcept = default;

        size_type size_bytes()const noexcept { return sizeof(element_type) * size(); }

        const iterator begin ()const noexcept { return iterator{ base_type::_Mydata }; }
        const iterator end   ()const noexcept { return iterator{ base_type::_Mydata + size() }; }
        const iterator rbegin()const noexcept { return reverse_iterator(end()); }
        const iterator rend  ()const noexcept { return reverse_iterator(begin()); }

        const reference front()const noexcept { return *base_type::_Mydata; }
        const reference back ()const noexcept { return *(base_type::_Mydata + size() - 1); }

        const reference operator[](size_type i)const noexcept { ASSERT(i >= 0 && i < size()); return base_type::_Mydata[i]; }

        const bool empty()const noexcept { return size() == 0; }

        template<size_type Count>
        const Span<T, Count> first()const
        {
            ASSERT(Count <= size());
            return Span<T, Count>(base_type::_Mydata, Count);
        }
        const Span<T, dynamic_extent> first(size_type count)const
        {
            ASSERT(count <= size());
            return Span<T, dynamic_extent>(base_type::_Mydata, count);
        }

        template<size_type Count>
        const Span<T, Count> last()const
        {
            ASSERT(Count <= size());
            return Span<T, Count>(base_type::_Mydata + (size() - Count), Count);
        }
        const Span<T, dynamic_extent> last(size_type count)const
        {
            ASSERT(count <= size());
            return Span<T, dynamic_extent>(base_type::_Mydata + (size() - count), count);
        }

        template<size_type Offset, size_type Count>
        const Span<T, Count> subspan()const
        {
            ASSERT(Offset + Count <= size());
            return Span<T, Count>(base_type::_Mydata + Offset, Count);
        }
        const Span<T, dynamic_extent> subspan(size_type offset, size_type count)const
        {
            ASSERT(offset + count <= size());
            return Span<T, dynamic_extent>(base_type::_Mydata + offset, count);
        }
    };
}
#endif

#ifdef __cpp_lib_byte
#include <cstddef>
#endif
namespace EasyLib {
    template<typename T, std::size_t N>
    auto as_bytes(Span<T> s)
    {
#ifdef __cpp_lib_byte
        using byte = std::byte;
#else
        using byte = unsigned char;
#endif
        if constexpr (N == dynamic_extent)
            return Span<const byte>(
                reinterpret_cast<const byte*>(s.data()),
                s.size_bytes()
            );
        else
            return Span<const byte, sizeof(T)* N>(
                reinterpret_cast<const byte*>(s.data()),
                sizeof(T) * N
            );
    }
    template<typename T, std::size_t N>
    auto as_writable_bytes(Span<T> s)
    {
        static_assert(!std::is_const_v<T>, "type is not writable!");
#ifdef __cpp_lib_byte
        using byte = std::byte;
#else
        using byte = unsigned char;
#endif
        if constexpr (N == dynamic_extent)
            return Span<byte>(
                reinterpret_cast<byte*>(s.data()),
                s.size_bytes()
            );
        else
            return Span<byte, sizeof(T)* N>(
                reinterpret_cast<byte*>(s.data()),
                sizeof(T) * N
            );
    }
}
