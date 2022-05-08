#pragma once
#include <vector>
#include <span>

#include "Index.hpp"

namespace EasyLib {

    class MeshConnectivity
    {
    public:
        using ivec = std::vector<int_l>;

        MeshConnectivity();
        MeshConnectivity(const MeshConnectivity&) = default;
        MeshConnectivity& operator=(const MeshConnectivity&) = default;

        void clear()
        {
            size_ = 0;
            ia_   = nullptr;
            ja_   = nullptr;

            _ia_.clear();
            _ja_.clear();
        }

        void reserve(int max_nrow, int max_ndata)
        {
            _ia_.reserve(max_nrow + 1);
            _ja_.reserve(max_ndata);
        }

        std::span<int> push_back(int n, const int* data)
        {

        }
        std::span<int> push_back(int n)
        {

        }
        std::span<int> push_back(std::span<const int> data)
        {

        }
        std::span<int> push_back(const std::initializer_list<int>& list)
        {

        }

        template<size_t N>
        std::span<int> push_back(const int(&list)[N])
        {
            static_assert(N > 0, "invalid array length");
            return this->push_back(N, list);
        }

        inline int nrow()const { return static_cast<int>(ia_.size() - 1); }
        inline int ndata()const { return static_cast<int>(ja_.size()); }

        inline bool empty()const { return ia_.size() == 1; }

        inline const ivec& ia()const { return ia_; }
        inline const ivec& ja()const { return ja_; }

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
        int_l size_{ 0 };
        const int_l* ia_{ nullptr };
        const int_l* ja_{ nullptr };

        ivec _ia_;
        ivec _ja_;
    };
}
