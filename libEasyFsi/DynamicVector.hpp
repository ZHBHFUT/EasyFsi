#pragma once
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <span>

#include "Index.hpp"
#include "DynamicArray.hpp"

namespace EasyLib {

    class Communicator;

    class DynamicVector : public DynamicArray<double, 1>
    {
    public:
        using DynamicArray<double, 1>::DynamicArray;

        inline void zero()
        {
            std::memset(data(), 0, sizeof(element_type) * numel());
        }

        inline double norm_sq()const
        {
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            for (auto p = data(); p != data() + numel(); ++p)ret += (*p) * (*p);
            return ret;
        }

        inline double norm()const
        {
            return std::sqrt(norm_sq());
        }

        inline auto& min()const { return *std::min_element(data(), data() + numel()); }
        inline auto& min()      { return *std::min_element(data(), data() + numel()); }
        inline auto& max()const { return *std::max_element(data(), data() + numel()); }
        inline auto& max()      { return *std::max_element(data(), data() + numel()); }

        inline double mean()const
        {
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            for (auto p = data(); p != data() + numel(); ++p)ret += *p;
            return ret / static_cast<double>(numel());
        }

        inline double dot(const DynamicVector& other)const
        {
            if (numel() != other.numel())throw std::runtime_error("size not agree!");
            if (empty())throw std::runtime_error("vector is empty!");
            double ret = 0;
            auto p = data();
            auto q = other.data();
            for (; p != data() + numel(); ++p, ++q)ret += (*p) * (*q);
            return ret;
        }
    };

    inline double norm   (const DynamicVector& v) { return v.norm(); }
    inline double norm_sq(const DynamicVector& v) { return v.norm_sq(); }
    inline double dot    (const DynamicVector& a, const DynamicVector& b) { return a.dot(b); }
}
