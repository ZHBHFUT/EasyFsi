#pragma once
#include <array>
#include <vector>
#include <span>

#include "Index.hpp"

namespace EasyLib {

    class MeshConnectivity
    {
    public:
        using value_type = int_l;
        using ivec = std::vector<int_l>;

        MeshConnectivity() = default;
        MeshConnectivity(const MeshConnectivity&) = default;
        MeshConnectivity& operator=(const MeshConnectivity&) = default;

        MeshConnectivity(MeshConnectivity&& mc)noexcept;
        MeshConnectivity& operator = (MeshConnectivity&& mc)noexcept;

        void clear()noexcept;

        void reserve(int_l max_nrow, int_l max_ndata);

        std::span<int_l> push_back(int_l n, const int_l* data = nullptr);
        std::span<int_l> push_back(std::span<const int_l> data);
        std::span<int_l> push_back(const std::initializer_list<int_l>& list);

        template<size_t N>
        std::span<int_l> push_back(const int_l(&list)[N])
        {
            static_assert(N > 0, "invalid array length");
            return this->push_back(N, list);
        }

        template<size_t N>
        std::span<int_l> push_back(const std::array<int_l, N>& list)
        {
            static_assert(N > 0, "invalid array length");
            return this->push_back(N, list.data());
        }

        inline int_l nrow ()const noexcept { return static_cast<int_l>(ia_.size() - 1); }
        inline int_l ndata()const noexcept { return static_cast<int_l>(ja_.size()); }

        inline bool empty()const noexcept { return ia_.size() == 1; }

        inline const ivec& ia()const noexcept { return ia_; }
        inline const ivec& ja()const noexcept { return ja_; }

        inline std::span<int> operator[](int row)
        {
            return std::span<int>{
                ja_.data() + ia_.at(row),
                    static_cast<size_t>(ia_[row + 1] - ia_[row])
            };
        };
        inline std::span<const int> operator[](int row)const
        {
            return std::span<const int>{
                ja_.data() + ia_.at(row),
                    static_cast<size_t>(ia_[row + 1] - ia_[row])
            };
        };

        inline int operator()(int i, int j)const
        {
            return ja_.at(ia_.at(i) + j);
        }

        static void flip(const MeshConnectivity& a2b, int nb, MeshConnectivity& b2a);

        friend class Boundary;
        friend class Interpolator;
        friend class Communicator;

    private:
        ivec ia_{ 0 };
        ivec ja_;
    };
}
