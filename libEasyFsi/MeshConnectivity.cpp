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
//! @file       MeshConnectivity.cpp
//!             The implement of MeshConnectivity class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------
#include "MeshConnectivity.hpp"

#include <stdexcept>
#include <algorithm>

#include "Assert.hpp"
#include "Span.hpp"

namespace EasyLib {

    MeshConnectivity::MeshConnectivity(MeshConnectivity&& mc)noexcept
        :ia_(std::move(mc.ia_)),
        ja_(std::move(mc.ja_))
    {}

    MeshConnectivity& MeshConnectivity::operator=(MeshConnectivity&& mc)noexcept
    {
        ia_ = std::move(mc.ia_);
        ja_ = std::move(mc.ja_);
        return *this;
    }

    void MeshConnectivity::clear()noexcept
    {
        ia_.clear(); ia_.push_back(0);
        ja_.clear();
    }

    void MeshConnectivity::reserve(value_type max_nrow, value_type max_ndata)
    {
        ia_.reserve(max_nrow + 1);
        ja_.reserve(max_ndata);
    }

    Span<int_l> MeshConnectivity::push_back(value_type n, const value_type* data/* = nullptr*/)
    {
        if (n < 0)n = 0;

        auto pos = ia_.back();
        ia_.push_back(pos + n);
        ja_.resize(ia_.back());
        if (n > 0 && data)std::copy(data, data + n, ja_.begin() + pos);
        return n > 0
            ? Span<int>(ja_.data() + pos, n)
            : Span<int>{};
    }
    Span<int_l> MeshConnectivity::push_back(Span<const value_type> data)
    {
        return push_back(static_cast<int_l>(data.size()), data.data());
    }
    Span<int_l> MeshConnectivity::push_back(const std::initializer_list<value_type>& list)
    {
        auto pos = ia_.back();
        ia_.push_back(pos + static_cast<int_l>(list.size()));
        ja_.resize(ia_.back());
        if (list.size() > 0)std::copy(list.begin(), list.end(), ja_.begin() + pos);
        return list.size() > 0
            ? Span<int>(ja_.data() + pos, list.size())
            : Span<int>{};
    }
    void MeshConnectivity::flip(const MeshConnectivity& a2b, value_type nb, MeshConnectivity& b2a)
    {
        b2a.clear();

        ivec count(nb);
        for (auto b : a2b.ja_) {
            ASSERT(b >= 0 && b < nb);
            if (b < 0 || b >= nb)throw std::runtime_error("invalid index!");
            ++count[b];
        }
        int nmax = *std::max_element(count.begin(), count.end());
        ivec faces(nmax);

        b2a.ia_.resize((size_t)nb + 1, 0);
        b2a.ja_.resize(a2b.ja_.size());
        for (int i = 0; i < nb; ++i)b2a.ia_[i + 1] = b2a.ia_[i] + count[i];

        std::fill(count.begin(), count.end(), 0);
        const auto na = a2b.nrow();
        for (int a = 0; a < na; ++a) {
            auto beg = a2b.ia_[a + 0];
            auto end = a2b.ia_[a + 1];
            for (int k = beg; k < end; ++k) {
                auto b = a2b.ja_[k];
                b2a.ja_[b2a.ia_[b] + count[b]] = a;
                ++count[b];
            }
        }
    }

}
