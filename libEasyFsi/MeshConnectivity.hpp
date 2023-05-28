#pragma once
#include <array>
#include <span>
#include <vector>

#include "Index.hpp"

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

        std::span<value_type> push_back(value_type n, const value_type* data = nullptr);
        std::span<value_type> push_back(std::span<const value_type> data);
        std::span<value_type> push_back(const std::initializer_list<value_type>& list);

        template<size_t N>
        std::span<value_type> push_back(const value_type(&list)[N])
        {
            static_assert(N > 0, "invalid array length");
            return this->push_back(N, list);
        }

        template<size_t N>
        std::span<value_type> push_back(const std::array<value_type, N>& list)
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

        inline std::span<value_type> operator[](int row)
        {
            return std::span<value_type>{
                ja_.data() + ia_.at(row),
                    static_cast<size_t>(ia_[row + 1] - ia_[row])
            };
        };
        inline std::span<const value_type> operator[](int row)const
        {
            return std::span<const value_type>{
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
