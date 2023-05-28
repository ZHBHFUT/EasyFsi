#pragma once
#include <vector>
#include <unordered_map>

#include "Index.hpp"

namespace EasyLib {

    class Boundary;
    class Communicator;

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
            return g2l_.contains(global_id);
        }

        inline std::pair<bool, int_l> find(int_g global_id)const
        {
            auto it = g2l_.find(global_id);
            return it != g2l_.end()
                ? std::make_pair(true, it->second)
                : std::make_pair(true, invalid_id);
        }

        inline constexpr bool empty()const noexcept { return l2g_.empty(); }

        inline constexpr int_l size()const noexcept { return static_cast<int_l>(l2g_.size()); }

        inline constexpr const int_g* data()const { return l2g_.data(); }

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
