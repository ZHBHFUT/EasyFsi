#include <stdexcept>
#include <algorithm>

#include "Assert.hpp"
#include "MeshConnectivity.hpp"

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
        ia_.clear();
        ja_.clear();
    }

    void MeshConnectivity::reserve(int_l max_nrow, int_l max_ndata)
    {
        ia_.reserve(max_nrow + 1);
        ja_.reserve(max_ndata);
    }

    std::span<int_l> MeshConnectivity::push_back(int_l n, const int_l* data/* = nullptr*/)
    {
        if (n < 0)n = 0;

        auto pos = ia_.back();
        ia_.push_back(pos + n);
        ja_.resize(ia_.back());
        if (n > 0 && data)std::copy(data, data + n, ja_.begin() + pos);
        return n > 0
            ? std::span<int>(ja_.data() + pos, n)
            : std::span<int>{};
    }
    std::span<int_l> MeshConnectivity::push_back(std::span<const int> data)
    {
        return push_back(static_cast<int_l>(data.size()), data.data());
    }
    std::span<int_l> MeshConnectivity::push_back(const std::initializer_list<int>& list)
    {
        auto pos = ia_.back();
        ia_.push_back(pos + static_cast<int_l>(list.size()));
        ja_.resize(ia_.back());
        if (list.size() > 0)std::copy(list.begin(), list.end(), ja_.begin() + pos);
        return list.size() > 0
            ? std::span<int>(ja_.data() + pos, list.size())
            : std::span<int>{};
    }
    void MeshConnectivity::flip(const MeshConnectivity& a2b, int nb, MeshConnectivity& b2a)
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
