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
//! @file       DynamicVector.hpp
//!             The definition of DynamicVector class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "Index.hpp"
#include "DynamicArray.hpp"

namespace EasyLib {

    class DynamicVector : public DynamicArray<double, 1>
    {
    public:
        using DynamicArray<double, 1>::DynamicArray;

        inline void zero()
        {
            std::memset(data(), 0, sizeof(element_type) * numel());
        }

        inline double norm_sq()const
        {
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            for (auto p = data(); p != data() + numel(); ++p)ret += (*p) * (*p);
            return ret;
        }

        inline double norm()const
        {
            return std::sqrt(norm_sq());
        }

        inline auto& min()const { return *std::min_element(data(), data() + numel()); }
        inline auto& min()      { return *std::min_element(data(), data() + numel()); }
        inline auto& max()const { return *std::max_element(data(), data() + numel()); }
        inline auto& max()      { return *std::max_element(data(), data() + numel()); }

        inline double mean()const
        {
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            for (auto p = data(); p != data() + numel(); ++p)ret += *p;
            return ret / static_cast<double>(numel());
        }

        inline double dot(const DynamicVector& other)const
        {
            if (numel() != other.numel())throw std::runtime_error("size not agree!");
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            auto p = data();
            auto q = other.data();
            for (; p != data() + numel(); ++p, ++q)ret += (*p) * (*q);
            return ret;
        }

        friend class Communicator;
    };

    inline double norm   (const DynamicVector& v) { return v.norm(); }
    inline double norm_sq(const DynamicVector& v) { return v.norm_sq(); }
    inline double dot    (const DynamicVector& a, const DynamicVector& b) { return a.dot(b); }
}
