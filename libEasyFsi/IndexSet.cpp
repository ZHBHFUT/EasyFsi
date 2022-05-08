#include "IndexSet.hpp"

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
    IndexSet::IndexSet(const lvec& _l2g, const imap _g2l)
        :l2g_(_l2g), g2l_(_g2l)
    {
        //TODO: check
    }
    IndexSet::IndexSet(lvec&& _l2g, imap&& _g2l)
        :l2g_(std::move(_l2g)),
        g2l_(std::move(_g2l))
    {
        //TODO: check
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

    void IndexSet::create(int_l count, const int_g* global_ids)
    {
        clear();

        count = std::max(int_l{ 0 }, count);

        // build g2l
        int_l ndup = 0;
        for (int_l l = 0; l < count; ++l) {
            if (!g2l_.contains(global_ids[l]))
                g2l_.emplace(global_ids[l], l);
            else
                + ndup;
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

        g2l_.emplace(global_id, static_cast<int_l>(g2l_.size()));
        l2g_.emplace_back(global_id);
    }

    void IndexSet::build_g2l_()
    {
        g2l_.clear();
        for (int_l l = 0; auto g : l2g_) {
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
}