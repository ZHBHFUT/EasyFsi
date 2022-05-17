#pragma once
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <vector>

#include "Index.hpp"

namespace EasyLib {

    class Communicator;

    class DynamicVector : public std::vector<double>
    {
    public:
        using base_type = std::vector<double>;

        using std::vector<double>::vector;
        using std::vector<double>::operator=;
        using std::vector<double>::reserve;
        using std::vector<double>::resize;
        using std::vector<double>::push_back;

        inline static constexpr int dimension = 1;

        inline constexpr int ndim()const noexcept { return 1; }

        inline void fill(double value)
        {
            std::fill(begin(), end(), value);
        }

        inline void copy_from(const DynamicVector& vec)
        {
            if (vec.size() != size())resize(vec.size()); // throw std::runtime_error("size not agree!");
            std::copy(vec.begin(), vec.end(), this->begin());
        }

        inline void swap(DynamicVector& vec)noexcept
        {
            std::swap(static_cast<base_type&>(*this), static_cast<base_type&>(vec));
        }

        inline double norm()const
        {
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            for (auto x : *this)ret += x * x;
            return std::sqrt(ret);
        }

        inline double norm_sq()const
        {
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            for (auto x : *this)ret += x * x;
            return ret;
        }

        inline auto& min()const { return *std::min_element(begin(), end()); }
        inline auto& min()      { return *std::min_element(begin(), end()); }
        inline auto& max()const { return *std::max_element(begin(), end()); }
        inline auto& max()      { return *std::max_element(begin(), end()); }

        inline double mean()const
        {
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            for (auto x : *this)ret += x;
            return ret / static_cast<double>(size());
        }

        inline double dot(const DynamicVector& other)const
        {
            if (size() != other.size())throw std::runtime_error("size not agree!");
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            auto p = begin();
            auto q = other.begin();
            for (; p != end(); ++p, ++q)ret += (*p) * (*q);
            return ret;
        }

        inline       auto& operator[](int_l i)      noexcept { return base_type::operator[](i); }
        inline const auto& operator[](int_l i)const noexcept { return base_type::operator[](i); }

        friend double norm(const DynamicVector& v)
        {
            return v.norm();
        }
        friend double norm_sq(const DynamicVector& v)
        {
            return v.norm_sq();
        }
        friend double dot(const DynamicVector& a, const DynamicVector& b)
        {
            return a.dot(b);
        }

        friend class Communicator;
    };

}
