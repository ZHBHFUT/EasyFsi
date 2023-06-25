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
//! @file       MeshConnectivity.hpp
//!             The definition MeshConnectivity class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <array>
#include <vector>

#include "Index.hpp"
#include "Span.hpp"

namespace EasyLib {

    class MeshConnectivity
    {
    public:
        using value_type = int_l;
        using ivec       = std::vector<int_l>;

        MeshConnectivity() = default;
        MeshConnectivity(const MeshConnectivity&) = default;
        MeshConnectivity& operator=(const MeshConnectivity&) = default;

        MeshConnectivity(MeshConnectivity&& mc)noexcept;
        MeshConnectivity& operator = (MeshConnectivity&& mc)noexcept;

        void clear()noexcept;

        void reserve(value_type max_nrow, value_type max_ndata);

        Span<value_type> push_back(value_type n, const value_type* data = nullptr);
        Span<value_type> push_back(Span<const value_type> data);
        Span<value_type> push_back(const std::initializer_list<value_type>& list);

        template<size_t N>
        Span<value_type> push_back(const value_type(&list)[N])
        {
            static_assert(N > 0, "invalid array length");
            return this->push_back(N, list);
        }

        template<size_t N>
        Span<value_type> push_back(const std::array<value_type, N>& list)
        {
            static_assert(N > 0, "invalid array length");
            return this->push_back(N, list.data());
        }

        inline value_type nrow ()const noexcept { return static_cast<value_type>(ia_.size() - 1); }
        inline value_type ndata()const noexcept { return static_cast<value_type>(ja_.size()); }
        inline value_type ndata(value_type row)const noexcept { return static_cast<value_type>(ia_.at(row + 1) - ia_.at(row)); }

        inline bool empty()const noexcept { return ia_.size() == 1; }

        inline const ivec& ia()const noexcept { return ia_; }
        inline const ivec& ja()const noexcept { return ja_; }

        inline Span<value_type> operator[](int row)
        {
            return Span<value_type>{
                ja_.data() + ia_.at(row),
                    static_cast<size_t>(ia_[row + 1] - ia_[row])
            };
        };
        inline Span<const value_type> operator[](int row)const
        {
            return Span<const value_type>{
                ja_.data() + ia_.at(row),
                    static_cast<size_t>(ia_[row + 1] - ia_[row])
            };
        };

        inline value_type operator()(value_type i, value_type j)const
        {
            return ja_.at(ia_.at(i) + j);
        }

        static void flip(const MeshConnectivity& a2b, value_type nb, MeshConnectivity& b2a);

        friend class Boundary;
        friend class Interpolator;
        friend class Communicator;

    private:
        ivec ia_{ 0 };
        ivec ja_;
    };
}
