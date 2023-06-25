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
//! @file       DynamicMatrix.hpp
//!             The definition of DynamicMatrix class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <algorithm>
#include <stdexcept>
#include <array>
#include <vector>

#include "Assert.hpp"
#include "LinAlgs.h"
#include "Span.hpp"
#include "DynamicVector.hpp"

namespace EasyLib {

    class DynamicMatrix : public DynamicArray<double, 2>
    {
    public:
        using DynamicArray<double, 2>::DynamicArray;

        inline constexpr bool is_square()const noexcept { return extent(0) == extent(1); }

        inline size_type nrow()const noexcept { return extent(0); }
        inline size_type ncol()const noexcept { return extent(1); }

        inline void zero()
        {
            std::memset(data(), 0, sizeof(element_type) * numel());
        }

        inline void identity()
        {
            zero();
            auto n = std::min(extent(0), extent(1));
            for (size_type i = 0; i < n; ++i)
                this->operator()(i * extent(1) + i) = 1;
        }

        //! @brief inverse this matrix.
        //! @return true=OK，false=Failed(matrix is singular)
        inline bool inverse()
        {
            if (!is_square())throw std::logic_error("DynamicMatrix::inverse, matrix is not square!");
            DynamicArray<int, 1> buffer(2 * extent(0), 0);
            return inverse(buffer);
        }
        inline bool inverse(DynamicArray<int, 1>& buffer)
        {
            if (!is_square())throw std::logic_error("DynamicMatrix::inverse, matrix is not square!");
            buffer.resize(2 * extent(0));
            int singular = 0;
            mat_inverse((int)extent(0), data(), buffer.data(), &singular);
            return singular != 0;
        }

        //! @brief compute y = A.x
        inline void apply(const DynamicVector& x, DynamicVector& y)const
        {
            if (ncol() != x.numel() || nrow() != y.numel())
                throw std::invalid_argument("DynamicMatrix::apply(), size not agree!");
            mat_apply_vec((int)nrow(), (int)ncol(), data(), x.data(), y.data());
        }
        inline void apply(const DynamicMatrix& x, DynamicMatrix& y)const
        {
            if (ncol() != x.nrow() ||
                nrow() != y.nrow() ||
                x.ncol() != y.ncol())
                throw std::invalid_argument("DynamicMatrix::apply(), size not agree!");
            mat_apply_mat((int)nrow(), (int)ncol(), (int)x.ncol(), data(), x.data(), y.data());
        }

        //! @brief compute y = A . x + y
        inline void apply_add(const DynamicVector& x, DynamicVector& y)const
        {
            if (ncol() != x.numel() || nrow() != y.numel())
                throw std::invalid_argument("DynamicMatrix::mat_apply_add_vec(), size not agree!");
            mat_apply_add_vec((int)nrow(), (int)ncol(), data(), x.data(), y.data());
        }
        inline void apply_add(const DynamicMatrix& x, DynamicMatrix& y)const
        {
            if (ncol() != x.nrow() ||
                nrow() != y.nrow() ||
                x.ncol() != y.ncol())
                throw std::invalid_argument("DynamicMatrix::apply(), size not agree!");
            mat_apply_add_mat((int)nrow(), (int)ncol(), (int)x.ncol(), data(), x.data(), y.data());
        }

        friend class Communicator;
    };

}
