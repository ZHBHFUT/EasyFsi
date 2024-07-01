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
//! @file       IndexSet.hpp
//!             The definition IndexSet class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include <vector>
#include <unordered_map>

#include "Index.hpp"

namespace EasyLib {

    class Boundary;
    class Communicator;

    //! @brief IndexSet class store the local and global index mapping.
    class IndexSet
    {
    public:
        using local_index_type  = int_l;
        using global_index_type = int_g;

        using imap = std::unordered_map<int_g, int_l>;
        using lvec = std::vector<int_g>;

        IndexSet() = default;
        IndexSet(const IndexSet&) = default;
        IndexSet& operator=(const IndexSet&) = default;

        explicit IndexSet(const lvec& _l2g);
        explicit IndexSet(lvec&& _l2g);
        explicit IndexSet(const imap& _g2l);
        explicit IndexSet(imap&& _g2l);
        IndexSet(int_l count, const int_g* global_ids);
        IndexSet(const lvec& _l2g, const imap& _g2l);
        IndexSet(lvec&& _l2g, imap&& _g2l);

        IndexSet(IndexSet&& is)noexcept;
        IndexSet& operator=(IndexSet&& is)noexcept;

        void create(const lvec& global_ids);
        void create(lvec&& global_ids);
        void create(const imap& _g2l);
        void create(imap&& _g2l);
        void create(const lvec& _l2g, const imap& _g2l);
        void create(lvec&& _l2g, imap&& _g2l);
        void create(int_l count, const int_g* global_ids);

        void clear()noexcept;

        void reserve(int_l max_size);

        int_l add(int_g global_id);

        inline bool contains(int_g global_id)const noexcept
        {
#if __cplusplus>=202002L
            return g2l_.contains(global_id);
#else
            return g2l_.find(global_id) != g2l_.end();
#endif
        }

        inline std::pair<bool, int_l> find(int_g global_id)const
        {
            auto it = g2l_.find(global_id);
            return it != g2l_.end()
                ? std::make_pair(true, it->second)
                : std::make_pair(true, invalid_id);
        }

#if __cplusplus>=202002L
        inline constexpr bool empty()const noexcept { return l2g_.empty(); }
        inline constexpr int_l size()const noexcept { return static_cast<int_l>(l2g_.size()); }
        inline constexpr const int_g* data()const { return l2g_.data(); }
#else
        inline bool empty()const noexcept { return l2g_.empty(); }
        inline int_l size()const noexcept { return static_cast<int_l>(l2g_.size()); }
        inline const int_g* data()const { return l2g_.data(); }
#endif

        inline int_g operator[](int_l i)const { return l2g_[i]; }

        inline int_g l2g(int_l i)const { return l2g_.at(i); }
        inline int_l g2l(int_g i)const { return g2l_.at(i); }

        friend class Boundary;
        friend class Communicator;

    private:
        void build_g2l_();
        void build_l2g_();
        void check_()const;

    private:
        lvec l2g_;
        imap g2l_;
    };
}
