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

//!-----------------------------------------------------------------------------
//! \file VecMatAlg.hpp
//! \brief Linear algebraic algorithms for fixed-length vector and matrix.
//! \author ZHANG Bing, zhangbing@hfut.edu.cn
//! \date 2024-04-18
//! \copyright (c)2024 All rights reserved.
//!-----------------------------------------------------------------------------

#include <cmath>       // sqrt, hypot
#include <cstring>     // memcpy
#include <type_traits> //
#include <algorithm>   // min, max
#include <numeric>     // iota
#include <array>
#ifdef _MSC_VER
#include <crtdbg.h>
#else
#include <cassert>
#ifdef _DEBUG
#define _ASSERT(x) assert(x)
#else
#define _ASSERT(x) ((void)0)
#endif
#endif
#if __cplusplus >= 202002L
#include <concepts> // std::floating_point
#define FLOAT_T std::floating_point
#define CHECK_T
#else
#define FLOAT_T typename
#define CHECK_T , FLOAT_T std::enable_if_t<std::is_floating_point<T>::value, bool> = true
#endif

namespace EasyLib {
    namespace VecMatAlg {
        template<FLOAT_T T, std::size_t N CHECK_T>                    using tvec = std::array<T, N>;
        template<FLOAT_T T, std::size_t M, std::size_t N = M CHECK_T> using tmat = std::array<std::array<T, N>, M>;

        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr tvec<T, N> zeros()
        {
            return tvec<T, N>{0};
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr tmat<T, M, N> zeros()
        {
            return tvec<T, M, N>{0};
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr tvec<T, N> ones()
        {
            tvec<T, N> a; a.fill(1); return a;
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr tmat<T, M, N> ones()
        {
            tmat<T, M, N> a; a.fill(1); return a;
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr tmat<T, M, N> eye()
        {
            constexpr auto R = M >= N ? M : N;
            tmat<T, M, N> a{ 0 };
            for (std::size_t i = 0; i < R; ++i)
                a[i][i] = 1;
            return a;
        }

        template<FLOAT_T T CHECK_T>
        constexpr T cross(const tvec<T, 2>& a, const tvec<T, 2>& b)
        {
            return a[0] * b[1] - a[1] * b[0];
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 3> cross(const tvec<T, 3>& a, const tvec<T, 3>& b)
        {
            return tvec<T, 3>{
                a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0]
            };
        }

        template<FLOAT_T T CHECK_T>
        constexpr T dot(const tvec<T, 1>& a, const tvec<T, 1>& b)
        {
            return a[0] * b[0];
        }
        template<FLOAT_T T CHECK_T>
        constexpr T dot(const tvec<T, 2>& a, const tvec<T, 2>& b)
        {
            return a[0] * b[0] + a[1] * b[1];
        }
        template<FLOAT_T T CHECK_T>
        constexpr T dot(const tvec<T, 3>& a, const tvec<T, 3>& b)
        {
            return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr T dot(const tvec<T, N>& a, const tvec<T, N>& b)
        {
            T r{ 0 };
            for (std::size_t j = 0; j < N; ++j)
                r += a[0][j] * b[j];
            return r;
        }
        template<FLOAT_T T, std::size_t M, std::size_t N, std::size_t L CHECK_T>
        constexpr tmat<T, M, L> dot(const tmat<T, M, N>& a, const tmat<T, N, L>& b)
        {
            tmat<T, M, L> r;
            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t j = 0; j < L; ++j) {
                    T rij{ 0 };
                    for (std::size_t k = 0; k < N; ++k)
                        rij += a[i][k] * b[k][j];
                    r[i][j] = rij;
                }
            }
            return r;
        }
        template<FLOAT_T T, std::size_t N, std::size_t L CHECK_T>
        constexpr tvec<T, L> dot(const tmat<T, 1, N>& a, const tmat<T, N, L>& b)
        {
            tvec<T, L> r;
            for (std::size_t i = 0; i < L; ++i) {
                T ri{ 0 };
                for (std::size_t j = 0; j < N; ++j)
                    ri += a[0][j] * b[j][i];
                r[i] = ri;
            }
            return r;
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr tvec<T, M> dot(const tmat<T, M, N>& a, const tvec<T, N>& b)
        {
            tvec<T, M> r;
            for (std::size_t i = 0; i < M; ++i) {
                T ri{ 0 };
                for (std::size_t j = 0; j < N; ++j)
                    ri += a[i][j] * b[j];
                r[i] = ri;
            }
            return r;
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr tvec<T, N> dot(const tvec<T, M>& a, const tmat<T, M, N>& b)
        {
            tvec<T, N> r;
            for (std::size_t i = 0; i < N; ++i) {
                T ri{ 0 };
                for (std::size_t j = 0; j < M; ++j)
                    ri += a[j] * b[j][i];
                r[i] = ri;
            }
            return r;
        }
        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr T dot(const tvec<T, M>& a, const tmat<T, M, 1>& b)
        {
            T r{ 0 };
            for (std::size_t j = 0; j < M; ++j)
                r += a[j] * b[j][0];
            return r;
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr T dot(const tmat<T, 1, N>& a, const tvec<T, N>& b)
        {
            T r{ 0 };
            for (std::size_t j = 0; j < N; ++j)
                r += a[0][j] * b[j];
            return r;
        }

        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 1> operator +(const tvec<T, 1>& a, const tvec<T, 1>& b)
        {
            return tvec<T, 1>{a[0] + b[0]};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 2> operator +(const tvec<T, 2>& a, const tvec<T, 2>& b)
        {
            return tvec<T, 2>{a[0] + b[0], a[1] + b[1]};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 3> operator +(const tvec<T, 3>& a, const tvec<T, 3>& b)
        {
            return tvec<T, 3>{a[0] + b[0], a[1] + b[1], a[2] + b[2]};
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr tvec<T, N> operator +(const tvec<T, N>& a, const tvec<T, N>& b)
        {
            tvec<T, N> r;
            for (std::size_t i = 0; i < N; ++i)r[i] = a[i] + b[i];
            return r;
        }

        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 1> operator -(const tvec<T, 1>& a, const tvec<T, 1>& b)
        {
            return tvec<T, 1>{a[0] - b[0]};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 2> operator -(const tvec<T, 2>& a, const tvec<T, 2>& b)
        {
            return tvec<T, 2>{a[0] - b[0], a[1] - b[1]};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 3> operator -(const tvec<T, 3>& a, const tvec<T, 3>& b)
        {
            return tvec<T, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr tvec<T, N> operator -(const tvec<T, N>& a, const tvec<T, N>& b)
        {
            tvec<T, N> r;
            for (std::size_t i = 0; i < N; ++i)r[i] = a[i] - b[i];
            return r;
        }

        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 1> operator -(const tvec<T, 1>& a)
        {
            return tvec<T, 1>{-a[0]};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 2> operator -(const tvec<T, 2>& a)
        {
            return tvec<T, 2>{-a[0], -a[1]};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 3> operator -(const tvec<T, 3>& a)
        {
            return tvec<T, 3>{-a[0], -a[1], -a[2]};
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr tvec<T, N> operator -(const tvec<T, N>& a)
        {
            tvec<T, N> r;
            for (std::size_t i = 0; i < N; ++i)r[i] = -a[i];
            return r;
        }

        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 1> operator *(const tvec<T, 1>& a, T b)
        {
            return tvec<T, 1>{a[0] * b};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 1> operator *(T b, const tvec<T, 1>& a)
        {
            return tvec<T, 1>{a[0] * b};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 2> operator *(const tvec<T, 2>& a, T b)
        {
            return tvec<T, 2>{a[0] * b, a[1] * b};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 2> operator *(T b, const tvec<T, 2>& a)
        {
            return tvec<T, 2>{a[0] * b, a[1] * b};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 3> operator *(const tvec<T, 3>& a, T b)
        {
            return tvec<T, 3>{a[0] * b, a[1] * b, a[2] * b};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 3> operator *(T b, const tvec<T, 3>& a)
        {
            return tvec<T, 3>{a[0] * b, a[1] * b, a[2] * b};
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr tvec<T, N> operator *(const tvec<T, N>& a, T b)
        {
            tvec<T, N> r;
            for (std::size_t i = 0; i < N; ++i)r[i] = a[i] * b;
            return r;
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr tvec<T, N> operator *(T b, const tvec<T, N>& a)
        {
            tvec<T, N> r;
            for (std::size_t i = 0; i < N; ++i)r[i] = a[i] * b;
            return r;
        }

        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 1> operator /(const tvec<T, 1>& a, T b)
        {
            _ASSERT(b != 0);
            return tvec<T, 1>{a[0] / b};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 2> operator /(const tvec<T, 2>& a, T b)
        {
            _ASSERT(b != 0);
            return tvec<T, 2>{a[0] / b, a[1] / b};
        }
        template<FLOAT_T T CHECK_T>
        constexpr tvec<T, 3> operator /(const tvec<T, 3>& a, T b)
        {
            _ASSERT(b != 0);
            return tvec<T, 3>{a[0] / b, a[1] / b, a[2] / b};
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr tvec<T, N> operator /(const tvec<T, N>& a, T b)
        {
            _ASSERT(b != 0);
            tvec<T, N> r;
            for (std::size_t i = 0; i < N; ++i)r[i] = a[i] / b;
            return r;
        }

        //! @brief compute: a += b
        template<FLOAT_T T CHECK_T>
        constexpr void self_add(tvec<T, 1>& a, T b)
        {
            std::get<0>(a) += b;
        }
        //! @brief compute: a += b
        template<FLOAT_T T CHECK_T>
        constexpr void self_add(tvec<T, 1>& a, const tvec<T, 1>& b)
        {
            std::get<0>(a) += std::get<0>(b);
        }
        //! @brief compute: a += b
        template<FLOAT_T T CHECK_T>
        constexpr void self_add(tvec<T, 2>& a, const tvec<T, 2>& b)
        {
            std::get<0>(a) += std::get<0>(b);
            std::get<1>(a) += std::get<1>(b);
        }
        //! @brief compute: a += b
        template<FLOAT_T T CHECK_T>
        constexpr void self_add(tvec<T, 3>& a, const tvec<T, 3>& b)
        {
            std::get<0>(a) += std::get<0>(b);
            std::get<1>(a) += std::get<1>(b);
            std::get<2>(a) += std::get<2>(b);
        }
        //! @brief compute: a += b
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr void self_add(tvec<T, N>& a, const tvec<T, N>& b)
        {
            for (std::size_t i = 0; i < N; ++i)a[i] += b[i];
        }

        //! @brief compute: a -= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_sub(tvec<T, 1>& a, T b)
        {
            std::get<0>(a) -= b;
        }
        //! @brief compute: a -= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_sub(tvec<T, 1>& a, const tvec<T, 1>& b)
        {
            std::get<0>(a) -= std::get<0>(b);
        }
        //! @brief compute: a -= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_sub(tvec<T, 2>& a, const tvec<T, 2>& b)
        {
            std::get<0>(a) -= std::get<0>(b);
            std::get<1>(a) -= std::get<1>(b);
        }
        //! @brief compute: a -= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_sub(tvec<T, 3>& a, const tvec<T, 3>& b)
        {
            std::get<0>(a) -= std::get<0>(b);
            std::get<1>(a) -= std::get<1>(b);
            std::get<2>(a) -= std::get<2>(b);
        }
        //! @brief compute: a -= b
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr void self_sub(tvec<T, N>& a, const tvec<T, N>& b)
        {
            for (std::size_t i = 0; i < N; ++i)a[i] -= b[i];
        }

        //! @brief compute: a *= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_mul(tvec<T, 1>& a, T b)
        {
            std::get<0>(a) *= b;
        }
        //! @brief compute: a *= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_mul(tvec<T, 2>& a, T b)
        {
            std::get<0>(a) *= b;
            std::get<1>(a) *= b;
        }
        //! @brief compute: a *= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_mul(tvec<T, 3>& a, T b)
        {
            std::get<0>(a) *= b;
            std::get<1>(a) *= b;
            std::get<2>(a) *= b;
        }
        //! @brief compute: a *= b
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr void self_mul(tvec<T, N>& a, T b)
        {
            for (std::size_t i = 0; i < N; ++i)a[i] *= b;
        }

        //! @brief compute: a /= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_div(tvec<T, 1>& a, T b)
        {
            std::get<0>(a) /= b;
        }
        //! @brief compute: a /= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_div(tvec<T, 2>& a, T b)
        {
            std::get<0>(a) /= b;
            std::get<1>(a) /= b;
        }
        //! @brief compute: a /= b
        template<FLOAT_T T CHECK_T>
        constexpr void self_div(tvec<T, 3>& a, T b)
        {
            std::get<0>(a) /= b;
            std::get<1>(a) /= b;
            std::get<2>(a) /= b;
        }
        //! @brief compute: a /= b
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr void self_div(tvec<T, N>& a, T b)
        {
            for (std::size_t i = 0; i < N; ++i)a[i] /= b;
        }

        template<FLOAT_T T CHECK_T>
        constexpr T norm_sq(const tvec<T, 1>& a)
        {
            return a[0] * a[0];
        }
        template<FLOAT_T T CHECK_T>
        constexpr T norm_sq(const tvec<T, 2>& a)
        {
            return a[0] * a[0] + a[1] * a[1];
        }
        template<FLOAT_T T CHECK_T>
        constexpr T norm_sq(const tvec<T, 3>& a)
        {
            return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr T norm_sq(const tvec<T, N>& a)
        {
            T r{ 0 };
            for (std::size_t i = 0; i < N; ++i)r += a[i] * a[i];
            return r;
        }

        template<FLOAT_T T CHECK_T>
        constexpr T norm(const tvec<T, 1>& a)
        {
            return std::fabs(a[0]);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T norm(const tvec<T, 2>& a)
        {
            return std::hypot(a[0], a[1]);// std::sqrt(a[0] * a[0] + a[1] * a[1]);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T norm(const tvec<T, 3>& a)
        {
#if __cplusplus < 201703L
            return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
#else
            return std::hypot(a[0], a[1], a[2]);
#endif
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr T norm(const tvec<T, N>& a)
        {
            return std::sqrt(norm_sq(a));
        }

        template<FLOAT_T T CHECK_T>
        constexpr T norm_infinity(const tvec<T, 1>& a)
        {
            return std::fabs(a[0]);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T norm_infinity(const tvec<T, 2>& a)
        {
            return std::max(std::fabs(a[0]), std::fabs(a[1]));
        }
        template<FLOAT_T T CHECK_T>
        constexpr T norm_infinity(const tvec<T, 3>& a)
        {
            return std::max(std::fabs(a[0]), std::max(std::fabs(a[1]), std::fabs(a[2])));
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr T norm_infinity(const tvec<T, N>& a)
        {
            T amax = 0;
            for (auto& x : a)amax = std::max(amax, std::fabs(x));
            return amax;
        }

        template<FLOAT_T T CHECK_T>
        constexpr T normalize(tvec<T, 1>& a)
        {
            a[0] = a[0] != 0 ? std::copysign(T{ 1 }, a) : 0;
            return std::fabs(a[0]);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T normalize(tvec<T, 2>& a)
        {
            auto l = norm(a);
            if (l > 0) { a[0] /= l; a[1] /= l; }
            return l;
        }
        template<FLOAT_T T CHECK_T>
        constexpr T normalize(tvec<T, 3>& a)
        {
            auto l = norm(a);
            if (l > 0) { a[0] /= l; a[1] /= l; a[2] /= l; }
            return l;
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr T normalize(tvec<T, N>& a)
        {
            auto r = norm(a);
            for (std::size_t i = 0; i < N; ++i)a[i] /= r;
            return r;
        }

        template<FLOAT_T T CHECK_T>
        constexpr T distance(const tvec<T, 1>& a, const tvec<T, 1>& b)
        {
            return std::fabs(b[0] - a[0]);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T distance(const tvec<T, 2>& a, const tvec<T, 2>& b)
        {
            return norm(b - a);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T distance(const tvec<T, 3>& a, const tvec<T, 3>& b)
        {
            return norm(b - a);
        }

        template<FLOAT_T T CHECK_T>
        constexpr T distance_sq(const tvec<T, 1>& a, const tvec<T, 1>& b)
        {
            return (b[0] - a[0]) * (b[0] - a[0]);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T distance_sq(const tvec<T, 2>& a, const tvec<T, 2>& b)
        {
            return norm_sq(b - a);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T distance_sq(const tvec<T, 3>& a, const tvec<T, 3>& b)
        {
            return norm_sq(b - a);
        }

        template<FLOAT_T T CHECK_T>
        constexpr T min(const tvec<T, 1>& a)
        {
            return std::get<0>(a);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T min(const tvec<T, 2>& a)
        {
            return std::min(std::get<0>(a), std::get<1>(a));
        }
        template<FLOAT_T T CHECK_T>
        constexpr T min(const tvec<T, 3>& a)
        {
            return std::min(std::get<0>(a), std::min(std::get<1>(a), std::get<2>(a)));
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr T min(const tvec<T, N>& a)
        {
            return *std::min_element(a.begin(), a.end());
        }
        template<FLOAT_T T CHECK_T>
        constexpr T max(const tvec<T, 1>& a)
        {
            return std::get<0>(a);
        }
        template<FLOAT_T T CHECK_T>
        constexpr T max(const tvec<T, 2>& a)
        {
            return std::max(std::get<0>(a), std::get<1>(a));
        }
        template<FLOAT_T T CHECK_T>
        constexpr T max(const tvec<T, 3>& a)
        {
            return std::max(std::get<0>(a), std::max(std::get<1>(a), std::get<2>(a)));
        }
        template<FLOAT_T T, std::size_t N CHECK_T>
        constexpr T max(const tvec<T, N>& a)
        {
            return *std::max_element(a.begin(), a.end());
        }

        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr void identity(tmat<T, M>& a)
        {
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t j = 0; j < M; ++j)
                    a[i][j] = i == j ? 1 : 0;
        }

        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr T det(const tmat<T, M>& a)
        {
            if constexpr (M == 1)
                return a[0][0];
            else if constexpr (M == 2)
                return a[0][0] * a[1][1] - a[0][1] * a[1][0];
            else if constexpr (M == 3) {
                auto b0 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
                auto b1 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
                auto b2 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
                return a[0][0] * b0 + a[0][1] * b1 + a[0][2] * b2;
            }
            else {
                tmat<T, M - 1> b;
                T d{ 0 }, s{ 1 };
                for (std::size_t col = 0; col < M; ++col) {
                    for (std::size_t i = 1; i < M; ++i) {
                        for (std::size_t j = 0, k = 0; j < M; ++j) {
                            if (j != col) {
                                b[i - 1][k] = a[i][j]; ++k;
                            }
                        }
                    }
                    d += s * a[0][col] * det(b);
                    s *= -1;
                }
                return d;
            }
        }

        template<FLOAT_T T CHECK_T>
        constexpr bool inverse(tmat<T, 2>& a)
        {
            auto det = a[0][0] * a[1][1] - a[0][1] * a[1][0]; _ASSERT(det != 0);
            auto t = a[0][0];
            a[0][0] = a[1][1] / det;
            a[1][1] = t / det;
            a[0][1] /= -det;
            a[1][0] /= -det;
            return det != 0;
        }
        template<FLOAT_T T CHECK_T>
        constexpr bool inverse(tmat<T, 3>& a)
        {
            T a0[3][3] = {
                {a[0][0],a[0][1],a[0][2]},
                {a[1][0],a[1][1],a[1][2]},
                {a[2][0],a[2][1],a[2][2]} };

            a[0][0] = a0[1][1] * a0[2][2] - a0[1][2] * a0[2][1];
            a[1][0] = a0[1][2] * a0[2][0] - a0[1][0] * a0[2][2];
            a[2][0] = a0[1][0] * a0[2][1] - a0[1][1] * a0[2][0];

            T det =
                  a0[0][0] * a[0][0]
                + a0[0][1] * a[1][0]
                + a0[0][2] * a[2][0]; _ASSERT(det != 0);

            a[0][1] = (a0[0][2] * a0[2][1] - a0[0][1] * a0[2][2]) / det;
            a[0][2] = (a0[0][1] * a0[1][2] - a0[0][2] * a0[1][1]) / det;
            a[1][1] = (a0[0][0] * a0[2][2] - a0[0][2] * a0[2][0]) / det;
            a[1][2] = (a0[0][2] * a0[1][0] - a0[0][0] * a0[1][2]) / det;
            a[2][1] = (a0[0][1] * a0[2][0] - a0[0][0] * a0[2][1]) / det;
            a[2][2] = (a0[0][0] * a0[1][1] - a0[0][1] * a0[1][0]) / det;
            a[0][0] /= det;
            a[1][0] /= det;
            a[2][0] /= det;
            return det != 0;
        }
        template<FLOAT_T T CHECK_T>
        constexpr bool inverse(tmat<T, 4>& a)
        {
            auto pa = &a[0][0];
            T m[16];
            std::memcpy(m, pa, sizeof(T) * 4 * 4);

            pa[ 0] =  m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
            pa[ 1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
            pa[ 2] =  m[1] * m[ 6] * m[15] - m[1] * m[ 7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[ 7] - m[13] * m[3] * m[6];
            pa[ 3] = -m[1] * m[ 6] * m[11] + m[1] * m[ 7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[ 9] * m[2] * m[ 7] + m[ 9] * m[3] * m[6];
            pa[ 4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
            pa[ 5] =  m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
            pa[ 6] = -m[0] * m[ 6] * m[15] + m[0] * m[ 7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[ 7] + m[12] * m[3] * m[6];
            pa[ 7] =  m[0] * m[ 6] * m[11] - m[0] * m[ 7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[ 8] * m[2] * m[ 7] - m[ 8] * m[3] * m[6];
            pa[ 8] =  m[4] * m[ 9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
            pa[ 9] = -m[0] * m[ 9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
            pa[10] =  m[0] * m[ 5] * m[15] - m[0] * m[ 7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[ 7] - m[12] * m[3] * m[5];
            pa[11] = -m[0] * m[ 5] * m[11] + m[0] * m[ 7] * m[ 9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[ 9] - m[ 8] * m[1] * m[ 7] + m[ 8] * m[3] * m[5];
            pa[12] = -m[4] * m[ 9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
            pa[13] =  m[0] * m[ 9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
            pa[14] = -m[0] * m[ 5] * m[14] + m[0] * m[ 6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[ 6] + m[12] * m[2] * m[5];
            pa[15] =  m[0] * m[ 5] * m[10] - m[0] * m[ 6] * m[ 9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[ 9] + m[ 8] * m[1] * m[ 6] - m[ 8] * m[2] * m[5];

            auto det = m[0] * pa[0] + m[1] * pa[4] + m[2] * pa[8] + m[3] * pa[12]; _ASSERT(det != 0);
            for (std::size_t i = 0; i < 16; ++i)pa[i] /= det;
            return det != 0;
        }
        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr bool inverse(tmat<T, M>& a)
        {
            std::array<std::size_t, M> is{ 0 }, js{ 0 };

            //
            // compute inversion of matrix using Gauss full pivot select method
            //

            for (std::size_t k = 0; k < M; ++k) {

                auto&& rk = a[k]; // ptr + k * M;// k-th row

                T d{ 0 };

                //pivot element
                for (std::size_t i = k; i < M; ++i) {
                    auto&& ri = a[i]; // a(i,:)
                    for (std::size_t j = k; j < M; j++) {
                        auto p = std::fabs(ri[j]);
                        if (p > d) {
                            d = p;
                            is[k] = i;
                            js[k] = j;
                        }
                    }
                }

                // matrix is singular if max pivot element is zero
                if (d == 0)return false;

                if (is[k] != k) {
                    auto&& rk2 = a[is[k]];// a(is[k],:)
                    for (std::size_t j = 0; j < M; ++j) {
                        std::swap(rk[j], rk2[j]);
                    }
                }

                if (js[k] != k) {
                    auto j = js[k];
                    for (std::size_t i = 0; i < M; ++i)
                        std::swap(a[i][k], a[i][j]);
                }

                rk[k] = T{ 1 } / rk[k];
                for (std::size_t j = 0; j < M; ++j) {
                    if (j != k)rk[j] *= rk[k];
                }

                for (std::size_t i = 0; i < M; ++i) {
                    if (i != k) {
                        auto&& ri = a[i]; // i-th row
                        for (std::size_t j = 0; j < M; ++j) {
                            if (j != k)ri[j] -= ri[k] * rk[j];
                        }
                    }
                }

                for (std::size_t i = 0; i < M; ++i) {
                    if (i != k)a[i][k] *= -rk[k];
                }
            }

            for (std::size_t k = M - 1; ; --k) {
                auto&& rk = a[k]; // a(k,:)
                if (js[k] != k) {
                    auto&& rk2 = a[js[k]]; // a(js[k],:);
                    for (std::size_t j = 0; j < M; j++) {
                        std::swap(rk[j], rk2[j]);
                    }
                }
                if (is[k] != k) {
                    auto k2 = is[k];
                    for (std::size_t i = 0; i < M; i++) {
                        std::swap(a[i][k], a[i][k2]);
                    }
                }
                if (k == 0)break;
            }

            return true;
        }

        template<FLOAT_T T CHECK_T>
        constexpr bool inverse_symm(tmat<T, 2>& a)
        {
            auto a00 = a[0][0];
            auto a01 = a[0][1];
            auto a11 = a[1][1];
            auto det = -a01 * a01 + a00 * a11; _ASSERT(det != 0);
            a[0][0] =  a11 / det;
            a[1][0] = a[0][1] = -a01 / det;
            a[1][1] =  a00 / det;
            return det != 0;
        }
        template<FLOAT_T T CHECK_T>
        constexpr bool inverse_symm(tmat<T, 3>& a)
        {
            auto a00 = a[0][0];
            auto a01 = a[0][1];
            auto a02 = a[0][2];
            auto a11 = a[1][1];
            auto a12 = a[1][2];
            auto a22 = a[2][2];

            auto det =
                2*a01 * a02 * a12
                - a02 * a02 * a11
                - a00 * a12 * a12
                - a01 * a01 * a22
                + a00 * a11 * a22; _ASSERT(det != 0);
            a[0][0] = (a11 * a22 - a12 * a12) / det;
            a[1][1] = (a00 * a22 - a02 * a02) / det;
            a[2][2] = (a00 * a11 - a01 * a01) / det;
            a[1][0] = a[0][1] = (a02 * a12 - a01 * a22) / det;
            a[2][0] = a[0][2] = (a01 * a12 - a02 * a11) / det;
            a[2][1] = a[1][2] = (a01 * a02 - a00 * a12) / det;
            return det != 0;
        }
        template<FLOAT_T T CHECK_T>
        constexpr bool inverse_symm(tmat<T, 4>& a)
        {
            auto a00 = a[0][0];
            auto a01 = a[0][1];
            auto a02 = a[0][2];
            auto a03 = a[0][3];
            auto a11 = a[1][1];
            auto a12 = a[1][2];
            auto a13 = a[1][3];
            auto a22 = a[2][2];
            auto a23 = a[2][3];
            auto a33 = a[3][3];

            a[0][0] =  a11 * a22 * a33 - a11 * a23 * a23 - a12 * a12 * a33 + a12 * a13 * a23 + a13 * a12 * a23 - a13 * a13 * a22;
            a[0][1] = -a01 * a22 * a33 + a01 * a23 * a23 + a12 * a02 * a33 - a12 * a03 * a23 - a13 * a02 * a23 + a13 * a03 * a22;
            a[0][2] =  a01 * a12 * a33 - a01 * a13 * a23 - a11 * a02 * a33 + a11 * a03 * a23 + a13 * a02 * a13 - a13 * a03 * a12;
            a[0][3] = -a01 * a12 * a23 + a01 * a13 * a22 + a11 * a02 * a23 - a11 * a03 * a22 - a12 * a02 * a13 + a12 * a03 * a12;
            a[1][1] =  a00 * a22 * a33 - a00 * a23 * a23 - a02 * a02 * a33 + a02 * a03 * a23 + a03 * a02 * a23 - a03 * a03 * a22;
            a[1][2] = -a00 * a12 * a33 + a00 * a13 * a23 + a01 * a02 * a33 - a01 * a03 * a23 - a03 * a02 * a13 + a03 * a03 * a12;
            a[1][3] =  a00 * a12 * a23 - a00 * a13 * a22 - a01 * a02 * a23 + a01 * a03 * a22 + a02 * a02 * a13 - a02 * a03 * a12;
            a[2][2] =  a00 * a11 * a33 - a00 * a13 * a13 - a01 * a01 * a33 + a01 * a03 * a13 + a03 * a01 * a13 - a03 * a03 * a11;
            a[2][3] = -a00 * a11 * a23 + a00 * a13 * a12 + a01 * a01 * a23 - a01 * a03 * a12 - a02 * a01 * a13 + a02 * a03 * a11;
            a[3][3] =  a00 * a11 * a22 - a00 * a12 * a12 - a01 * a01 * a22 + a01 * a02 * a12 + a02 * a01 * a12 - a02 * a02 * a11;

            auto det = a00 * a[0][0] + a01 * a[0][1] + a02 * a[0][2] + a03 * a[0][3]; _ASSERT(det != 0);
            a[0][0] /= det;
            a[0][1] /= det;
            a[0][2] /= det;
            a[0][3] /= det;
            a[1][1] /= det;
            a[1][2] /= det;
            a[1][3] /= det;
            a[2][2] /= det;
            a[2][3] /= det;
            a[3][3] /= det;

            a[1][0] = a[0][1];
            a[2][0] = a[0][2];
            a[2][1] = a[1][2];
            a[3][0] = a[0][3];
            a[3][1] = a[1][3];
            a[3][2] = a[2][3];

            return det != 0;
        }

        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr tmat<T, N, M> transpose(const tmat<T, M, N>& a)
        {
            tmat<T, N, M> aT;
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t j = 0; j < N; ++j)
                    aT[j][i] = a[i][j];
            return aT;
        }

        //! @brief Compute C = A^T . B
        template<FLOAT_T T, std::size_t M, std::size_t N, std::size_t L CHECK_T>
        constexpr tmat<T, N, L> transpose_dot(const tmat<T, M, N>& a, const tmat<T, M, L>& b)
        {
            tmat<T, N, L> c;
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < L; ++j) {
                    c[i][j] = T{ 0 };
                    for (std::size_t k = 0; k < M; ++k)
                        c[i][j] += a[k][i] * b[k][j];
                }
            return c;
        }

        //! @brief Compute C = A^T . B
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr tvec<T, N> transpose_dot(const tmat<T, M, N>& a, const tvec<T, M>& b)
        {
            tvec<T, N> c;
            for (std::size_t i = 0; i < N; ++i) {
                c[i] = T{ 0 };
                for (std::size_t k = 0; k < M; ++k)
                    c[i] += a[k][i] * b[k]; // a(:,i)*b(:)
            }
            return c;
        }

        //! @brief Do pivoted LU decomposition of a non-singular matrix.
        //! @tparam T element type of the matrix, should be floating number.
        //! @tparam M rank of the matrix.
        //! @param a     Input matrix, output as the results.
        //! @param perm  The permutation of row.
        //! @return true=succeeded, false=matrix is singular.
        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr bool lu_dcmp(tmat<T, M>& a, std::array<std::size_t, M>& perm/*, int& flag*/)
        {
            constexpr auto rank = M;

            T    absA{ 0 };
            bool succeed = true;

            //flag = 1; // No row interchanges yet.

            // Unit permutation matrix, P[N] initialized with N
            std::iota(perm.begin(), perm.end(), std::size_t{ 0 });
            
            for (std::size_t i = 0; i < rank; i++) {
                std::size_t imax = i;
                T maxA{ 0 };

                // Search largest pivot element of i-column
                for (std::size_t k = i; k < rank; k++)
                    if ((absA = std::fabs(a[k][i])) > maxA) {
                        maxA = absA;
                        imax = k;
                    }

                // failure, matrix is degenerate if largest pivot element is zero
                if (maxA == 0)succeed = false;

                if (imax != i) {
                    // pivoting P
                    std::swap(perm[i], perm[imax]);

                    // pivoting rows of A
                    for (std::size_t k = 0; k < rank; ++k)
                        std::swap(a[imax][k], a[i][k]);

                    // counting pivots starting from N (for determinant)
                    //flag = -flag;
                }

                auto aii = a[i][i];
                for (std::size_t j = i + 1; j < rank; j++) {
                    a[j][i] /= aii;

                    for (std::size_t k = i + 1; k < rank; k++)
                        a[j][k] -= a[j][i] * a[i][k];
                }
            }

            _ASSERT(succeed);
            return succeed;
        }

        //! @brief Solve linear systems using LU decomposition of A: LU.x = b
        //! @param [in]  lu 
        //! @param [in]  perm 
        //! @param [in]  b 
        //! @param [out] x 
        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr void lu_backsubs(const tmat<T, M>& lu, const std::array<std::size_t, M>& perm, const tvec<T, M>& b, tvec<T, M>& x)
        {
            constexpr auto rank = M;

            // back substitution: L * x = b
            for (std::size_t i = 0; i < rank; ++i) {
                //check perm
                //if (perm[i] < 0 || perm[i] >= rank)
                //    throw std::runtime_error("***ERROR*** lu_backsubs(), perm index out of range.");
                _ASSERT(perm[i] >= 0 && perm[i] < rank);

                x[i] = b[perm[i]];
                for (std::size_t j = 0; j < i; ++j)
                    x[i] -= lu[i][j] * x[j];
            }
            // back substitution: U * x = x
            for (std::size_t i = rank - 1; ; --i) {
                for (std::size_t j = i + 1; j < rank; ++j)
                    x[i] -= lu[i][j] * x[j];
                x[i] /= lu[i][i];
                if (i == 0)break;
            }
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr void lu_backsubs(const tmat<T, M>& lu, const std::array<std::size_t, M>& perm, const tmat<T, M, N>& b, tmat<T, M, N>& x)
        {
            constexpr auto rank = M;

            //check perm
            //for (std::size_t i = 0; i < rank; ++i)
            //    if (perm[i] < 0 || perm[i] >= rank)
            //        throw std::runtime_error("***ERROR*** lu_backsubs(), perm index out of range.");

            // Do back substitution.
            // 
            // loop each column of X
            for (std::size_t j = 0; j < N; ++j) {

                // back substitution: L * x = b ==> x(:,j)
                for (std::size_t i = 0; i < rank; ++i) {
                    _ASSERT(perm[i] >= 0 && perm[i] < rank);
                    T xi = b[perm[i]][j];
                    for (std::size_t k = 0; k < i; k++)
                        xi -= lu[i][k] * x[k][j];//xi-=lu(i,k)*x(l,j)
                    x[i][j] = xi;
                }

                // back substitution: U * x = x ==> x(:,j)
                for (std::size_t i = rank - 1; ; --i) {
                    T xij = x[i][j];// x(i,j)
                    for (std::size_t k = i + 1; k < rank; k++)
                        xij -= lu[i][k] * x[k][j]; // A(i, k) * x(k, j)
                    xij /= lu[i][i]; // A(i, i);
                    x[i][j] = xij;
                    if (i == 0)break;
                }
            }
        }

        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr bool lu_solve(const tmat<T, M>& a, const tvec<T, M>& b, tvec<T, M>& x)
        {
            std::array<std::size_t, M> perm{ 0 };
            tmat<T, M> lu = a;
            int flag = 0;

            if (!lu_dcmp(lu, perm, flag))return false;
            lu_backsubs(lu, perm, b, x);
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr bool lu_solve(const tmat<T, M>& a, const tmat<T, M, N>& b, tmat<T, M, N>& x)
        {
            std::array<std::size_t, M> perm{ 0 };
            tmat<T, M> lu = a;
            int flag = 0;

            if (!lu_dcmp(lu, perm, flag))return false;
            lu_backsubs(lu, perm, b, x);
        }

        //! @brief Compute the Cholesky decomposition of a symmetric positive matrix.
        //! @tparam T element type of the matrix, should be floating number.
        //! @tparam M rank of the matrix.
        //! @param [in]  a  The symmetric positive matrix need decomposition.
        //! @param [out] L  The lower-part of the results matrix. i.e. L.L' = a
        //! @return true=succeeded, false=matrix is not positive.
        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr bool chol_dcmp(const tmat<T, M>& a, tmat<T, M>& L)
        {
            constexpr auto rank = M;
            // compute lower part
            for (std::size_t i = 0; i < rank; ++i) {
                for (std::size_t j = 0; j <= i; ++j) {
                    T sum{ 0 };

                    // diagonal term
                    if (j == i) {
                        for (int k = 0; k < j; ++k)
                            sum += L[i][k] * L[i][k];
                        if (a[i][i] <= sum) {
                            _ASSERT(false);
                            return false;
                        }
                        L[i][i] = std::sqrt(a[i][i] - sum);
                    }
                    // off-diag of lower part
                    else {
                        for (std::size_t k = 0; k < j; ++k)
                            sum += L[i][k] * L[j][k];
                        L[i][j] = (a[i][j] - sum) / L[j][j];
                    }
                }
            }
            // set upper part to ZERO
            for (std::size_t i = 0; i < rank; ++i) {
                for (std::size_t j = i + 1; j < rank; ++j) {
                    L[i][j] = 0;
                }
            }
            return true;
        }

        //! @brief Compute the Cholesky decomposition of a symmetric positive matrix.
        //! @tparam T element type of the matrix, should be floating number.
        //! @tparam M rank of the matrix.
        //! @param [inout] a Input as the symmetric positive matrix, output as the lower part of the decomposition.
        //! @return true=succeeded, false=failed.
        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr bool chol_dcmp(tmat<T, M>& a)
        {
            tmat<T, M> L;
            auto r = chol_dcmp(a, L);
            a = L;
            return r;
        }

        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr void chol_backsubs(tmat<T, M>& L, const tvec<T, M>& b, tvec<T, M>& x)
        {
            constexpr auto rank = M;

            // L L' x = b
            // step-1: solve  L y = b
            tvec<T, M> y;
            for (std::size_t i = 0; i < rank; ++i) {
                T sum{ 0 };
                for (std::size_t j = 0; j < i; ++j)
                    sum += L[i][j] * y[j];
                y[i] = (b[i] - sum) / L[i][i];
            }
            // step-2: solve L' x = y
            for (std::size_t i = rank - 1;; --i) {
                T sum{ 0 };
                for (std::size_t j = i + 1; j < rank; ++j)
                    sum += L[j][i] * x[j];
                x[i] = (y[i] - sum) / L[i][i];
                if (i == 0)break;
            }
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr void chol_backsubs(tmat<T, M>& L, const tmat<T, M, N>& b, tmat<T, M, N>& x)
        {
            constexpr auto rank = M;

            // L L' x = b
            tvec<T, M> y;
            for (std::size_t k = 0; k < N; ++k) {
                // step-1: solve  L y = b
                for (std::size_t i = 0; i < rank; ++i) {
                    T sum{ 0 };
                    for (std::size_t j = 0; j < i; ++j)
                        sum += L[i][j] * y[j];
                    y[i] = (b[i][k] - sum) / L[i][i];
                }
                // step-2: solve L' x = y
                for (std::size_t i = rank - 1;; --i) {
                    T sum{ 0 };
                    for (std::size_t j = i + 1; j < rank; ++j)
                        sum += L[j][i] * x[j][k];
                    x[i][k] = (y[i] - sum) / L[i][i];
                    if (i == 0)break;
                }
            }
        }

        template<FLOAT_T T, std::size_t M CHECK_T>
        constexpr bool chol_solve(const tmat<T, M>& a, const tvec<T, M>& b, tvec<T, M>& x)
        {
            tmat<T, M> L;
            if (!chol_dcmp(a, L))return false;
            chol_backsubs(L, b, x);
            return true;
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        constexpr bool chol_solve(const tmat<T, M>& a, const tmat<T, M, N>& b, tmat<T, M, N>& x)
        {
            tmat<T, M> L = a;
            if (!chol_dcmp(L))return false;
            chol_backsubs(L, b, x);
            return true;
        }

        template<FLOAT_T T CHECK_T>
        tvec<T, 1>& view1(T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<tvec<T, 1>*>(ptr + offset);
        }
        template<FLOAT_T T CHECK_T>
        const tvec<T, 1>& view1(const T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<const tvec<T, 1>*>(ptr + offset);
        }
        template<FLOAT_T T CHECK_T>
        tvec<T, 2>& view2(T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<tvec<T, 2>*>(ptr + offset);
        }
        template<FLOAT_T T CHECK_T>
        const tvec<T, 2>& view2(const T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<const tvec<T, 2>*>(ptr + offset);
        }
        template<FLOAT_T T CHECK_T>
        tvec<T, 3>& view3(T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<tvec<T, 3>*>(ptr + offset);
        }
        template<FLOAT_T T CHECK_T>
        const tvec<T, 3>& view3(const T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<const tvec<T, 3>*>(ptr + offset);
        }

        template<FLOAT_T T, std::size_t ND CHECK_T>
        tvec<T, ND>& vec_view(T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<tvec<T, ND>*>(ptr + offset);
        }
        template<FLOAT_T T, std::size_t ND CHECK_T>
        const tvec<T, ND>& vec_view(const T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<const tvec<T, ND>*>(ptr + offset);
        }

        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        tmat<T, M, N>& mat_view(T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<tmat<T, M, N>*>(ptr + offset);
        }
        template<FLOAT_T T, std::size_t M, std::size_t N CHECK_T>
        const tmat<T, M, N>& mat_view(const T* ptr, std::size_t offset = 0)
        {
            return *reinterpret_cast<const tmat<T, M, N>*>(ptr + offset);
        }
    }
}
