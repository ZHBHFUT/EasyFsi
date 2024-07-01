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
//! @file       IndexSet.cpp
//!             The implement of IndexSet class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include "IndexSet.hpp"
#include <stdexcept>
#include <algorithm>

namespace EasyLib {

    IndexSet::IndexSet(const lvec& _l2g)
        :l2g_(_l2g)
    {
        build_g2l_();
        if (g2l_.size() != l2g_.size())build_l2g_();
    }
    IndexSet::IndexSet(lvec&& _l2g)
        :l2g_(std::move(_l2g))
    {
        build_g2l_();
        if (g2l_.size() != l2g_.size())build_l2g_();
    }
    IndexSet::IndexSet(const imap& _g2l)
        :g2l_(_g2l)
    {
        build_l2g_();
    }
    IndexSet::IndexSet(imap&& _g2l)
        :g2l_(std::move(_g2l))
    {
        build_l2g_();
    }
    IndexSet::IndexSet(int_l count, const int_g* global_ids)
    {
        create(count, global_ids);
    }
    IndexSet::IndexSet(const lvec& _l2g, const imap& _g2l)
        :l2g_(_l2g), g2l_(_g2l)
    {
        check_();
    }
    IndexSet::IndexSet(lvec&& _l2g, imap&& _g2l)
        :l2g_(std::move(_l2g)),
        g2l_(std::move(_g2l))
    {
        check_();
    }

    IndexSet::IndexSet(IndexSet&& is)noexcept
        :l2g_(std::move(is.l2g_)),
        g2l_(std::move(is.g2l_))
    {}

    IndexSet& IndexSet::operator=(IndexSet&& is)noexcept
    {
        l2g_ = std::move(is.l2g_);
        g2l_ = std::move(is.g2l_);
        return *this;
    }

    void IndexSet::create(const lvec& global_ids)
    {
        clear();
        create(static_cast<int_l>(global_ids.size()), global_ids.data());
    }
    void IndexSet::create(lvec&& global_ids)
    {
        clear();
        l2g_ = std::move(global_ids);
        build_g2l_();
        if (g2l_.size() != l2g_.size())build_l2g_();
    }
    void IndexSet::create(const imap& _g2l)
    {
        clear();
        g2l_ = _g2l;
        build_l2g_();
    }
    void IndexSet::create(imap&& _g2l)
    {
        clear();
        g2l_ = std::move(_g2l);
        build_l2g_();
    }
    void IndexSet::create(const lvec& _l2g, const imap& _g2l)
    {
        l2g_ = _l2g;
        g2l_ = _g2l;
        check_();
    }
    void IndexSet::create(lvec&& _l2g, imap&& _g2l)
    {
        l2g_ = std::move(_l2g);
        g2l_ = std::move(_g2l);
        check_();
    }
    void IndexSet::create(int_l count, const int_g* global_ids)
    {
        clear();

        count = std::max(int_l{ 0 }, count);

        // build g2l
        int_l ndup = 0;
        for (int_l l = 0; l < count; ++l) {
            if (!this->contains(global_ids[l]))
                g2l_.emplace(global_ids[l], l);
            else
                ++ndup;
        }

        // no duplicated ids
        if (!ndup) {
            l2g_.resize(count);
            std::copy(global_ids, global_ids + count, l2g_.data());
        }
        // duplicated id
        else {
            l2g_.resize(g2l_.size());
            for (auto& p : g2l_)
                l2g_[p.second] = p.first;
            //throw std::runtime_error("duplicated index detected!");
        }
    }

    void IndexSet::clear()noexcept
    {
        g2l_.clear();
        l2g_.clear();
    }

    void IndexSet::reserve(int_l max_size)
    {
        l2g_.reserve(max_size);
    }

    int_l IndexSet::add(int_g global_id)
    {
        auto it = g2l_.find(global_id);
        if (it != g2l_.end())return it->second;

        auto id = static_cast<int_l>(g2l_.size());
        g2l_.emplace(global_id, id);
        l2g_.emplace_back(global_id);
        return id;
    }

    void IndexSet::build_g2l_()
    {
        g2l_.clear();
        int_l l = 0;
        for (auto g : l2g_) {
            g2l_.try_emplace(g, l);
            ++l;
        }
    }
    void IndexSet::build_l2g_()
    {
        l2g_.resize(g2l_.size());
        for (auto& p : g2l_)
            l2g_[p.second] = p.first;
    }
    void IndexSet::check_()const
    {
        if (g2l_.size() != l2g_.size())
            throw std::runtime_error("size not agree!");

        for (const auto& p : g2l_) {
            if (p.second < 0 || p.second >= l2g_.size())
                throw std::runtime_error("local index out of range!");

            if (l2g_[p.second] != p.first)
                throw std::runtime_error("index not agree!");
        }
    }
}