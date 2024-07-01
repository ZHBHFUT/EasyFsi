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
//! \file InverseIsoMapping.hpp
//! \brief Inverse-Isoparametric-mapping algorithm for standard 1-D/2-D/3-D finite elements.
//! \author ZHANG Bing, zhangbing@hfut.edu.cn
//! \date 2024-04-18
//! \copyright (c)2024 All rights reserved.
//!-----------------------------------------------------------------------------
#pragma once
#include <limits>
#include <algorithm>

#include "VecMatAlg.hpp"
#include "ShapeFunctions.hpp"
#include "LevenbergMarquardt.hpp"

namespace EasyLib {

    //-----------------------------------------------------------------------
    //    Compute inversion isoparametric mappings for standard elements
    //-----------------------------------------------------------------------
    // 
    //  X(s,t) = {N1(s,t),N2(s,t),...}.{{xI,yI},{xJ,yJ},...}
    //         =          N({s,t})    .       Xe
    //  for given p={x,y}, find {s,t}:
    //      X({s,t}) - p = 0
    // 
    // ==>  f({s,t}) = X({s,t}) - p
    // 
    // using Newton method:
    //      {s,t}_{k+1} = {s,t}_k - J^-1 . f({s,t}_k)
    // where,
    //          {dNds} 
    //      J = {    } . Xe   !!! NOTE: J maybe singular for many cases !!!
    //          {dNdt}
    // 
    // Because the Jacobian matrix tends to become singular in many instances.
    // We employ optimization algorithms to solve this problem:
    //   Find a nearest location on the element for the point.
    //       min(dot(dot(Fn,Xe)-p, dot(Fn,Xe)-p)
    // 
    //-----------------------------------------------------------------------
    
    //! @brief Inverse isoparametric mappings for 1D linear bar element.
    //! @param Xe     Vertex array for BAR2: {xI,xJ}
    //! @param point  Query point: {x}
    //! @param t      Local coordinate
    //! @code
    //! double Xe[2] = {-1,1};
    //! double p = 0.2;
    //! double t;
    //! inverse_iso_mapping_bar2_1d(Xe,p,t);
    //! // ==> t = 0.2
    //! double q = 1/2*(1-t)*Xe[0] + 1/2*(1+t)*Xe[1]
    //! // ==>   = 0.2
    //! @endcode
    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_bar2_1d(const MatXe& Xe, const VecP& point, Real& t)
    {
        static constexpr auto eps = std::numeric_limits<Real>::epsilon();

        // We use analytic formula:
        //  +-----o-------+   -->  +-----o-------+
        //  xi    p       xj      -1             1

        auto d1 = Xe[1] - Xe[0];
        auto d2 = 2 * point - Xe[0] - Xe[1];
        _ASSERT(d1 != 0);
        t = d1 ? d2 / d1 : (d2 ? Real{ -1 } : Real{ -2 }); // d1*d2/(d1*d1)
    }
    //! @brief Inverse isoparametric mappings for 2D linear bar element.
    //! @param Xe     Vertex array for BAR2: {xI,xJ}
    //! @param point  Query point: {x,y}
    //! @param t      Local coordinate
    //! @code
    //! double Xe[2][2] = {-1,-1,1,1};
    //! double p[2] = {0.2,0.2};
    //! double t;
    //! inverse_iso_mapping_bar2_2d(Xe,p,t);
    //! // ==> t = 0.2
    //! double q[2] = {
    //!     1/2*(1-t)*Xe[0][0] + 1/2*(1+t)*Xe[1][0],
    //!     1/2*(1-t)*Xe[0][1] + 1/2*(1+t)*Xe[1][1]
    //! };
    //! // ==>   = {0.2,0.2}
    //! @endcode
    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_bar2_2d(const MatXe& Xe, const VecP& point, Real& t)
    {
        auto d1x = Xe[1][0] - Xe[0][0];
        auto d1y = Xe[1][1] - Xe[0][1];
        auto d2x = 2 * point[0] - Xe[1][0] - Xe[0][0];
        auto d2y = 2 * point[1] - Xe[1][1] - Xe[0][1];
        auto b = d1x * d1x + d1y * d1y;
        auto a = d1x * d2x + d1y * d2y;
        _ASSERT(b != 0);
        t = b ? a / b : (a ? Real{ -1 } : Real{ -2 });
    }
    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_bar2_3d(const MatXe& Xe, const VecP& point, Real& t)
    {
        auto d1x = Xe[1][0] - Xe[0][0];
        auto d1y = Xe[1][1] - Xe[0][1];
        auto d1z = Xe[1][2] - Xe[0][2];
        auto d2x = 2 * point[0] - Xe[1][0] - Xe[0][0];
        auto d2y = 2 * point[1] - Xe[1][1] - Xe[0][1];
        auto d2z = 2 * point[2] - Xe[1][2] - Xe[0][2];
        auto b = d1x * d1x + d1y * d1y + d1z * d1z;
        auto a = d1x * d2x + d1y * d2y + d1z * d2z;
        _ASSERT(b != 0);
        t = b ? a / b : (a ? Real{ -1 } : Real{ -2 });
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_bar3_1d(const MatXe& Xe, const VecP& point, Real& t, const int ix = 0)
    {
        // use analytic formula:
        // f(t) = x(t) - px
        //      = 1/2*(t * t - t)*xI + 1/2*(t * t + t)*xJ + (1 - t * t)*xK - px
        //      = (xI/2 + xJ/2 -xK)*t^2 + (xJ/2-xI/2)*t + xK-px
        //      = 0
        static constexpr auto eps = std::numeric_limits<Real>::epsilon();
        const auto a = Real{ 0.5 } * (Xe[0][ix] + Xe[1][ix]) - Xe[2][ix];
        const auto b = Real{ 0.5 } * (Xe[1][ix] - Xe[0][ix]);
        const auto c = Xe[2][ix] - point[ix];
        if     (std::fabs(a) >= eps) {
            const auto d = std::sqrt(b * b - 4 * a * c);
            auto t1 = -b + d / (2 * a);
            auto t2 = -b - d / (2 * a);
            t = t1 >= -1 && t1 <= 1 ? t1 : t2;
        }
        else if(std::fabs(b) >= eps)
            t = -c / b;
        else
            t = -1;

        _ASSERT(t >= -1 && t <= 1);
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_bar3_2d(const MatXe& Xe, const VecP& point, Real& t)
    {
        using namespace VecMatAlg;
        static constexpr int NN = 3;
        static constexpr int ND = 2;
        static constexpr int NC = 1;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;

        t = 0;
        MinimizeLM<NC, ND>::solve(
            [&Xe](const type x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_bar3(x, fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const type x, tvec<type, ND>& J) {
                tvec<type, NN> dfn;
                shape_function_diff_bar3(x, dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            point,
            t
        );
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_bar3_3d(const MatXe& Xe, const VecP& point, Real& t)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>; // std::remove_cvref_t<decltype(Xe[0][0])>
        static constexpr int NN = 3;
        static constexpr int ND = 3;
        static constexpr int NC = 1;

        t = 0;
        MinimizeLM<NC, ND>::solve(
            [&Xe](const type x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_bar3(x, fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const type x, tvec<type, ND>& J) {
                tvec<type, NN> dfn;
                shape_function_diff_bar3(x, dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            point,
            t
        );
    }
    
    //! @brief Inverse isoparametric mappings for triangle
    //! @param Xe Vertex array for TRI: {I,J,K}
    //! @param p  Query point: {x,y,z}
    //! @param LI Local coordinate for I
    //! @param LJ Local coordinate for J
    //! @param ix id for component-x, default is 0
    //! @param iy id for component-y, default is 1
    //! @code
    //! double Xe[4][0] = {
    //!     {0,0},
    //!     {1,0},
    //!     {0,1}};
    //! double p[2] = {0.2,0.3};
    //! double L1,L2;
    //! inverse_iso_mapping_tri3(Xe,p,L1,L2);
    //! // ==> L1=0.5, L2=0.2
    //! double L3=1-L1-L2; // ==> L3 = 0.3
    //! double q[2]={
    //!     L1 * Xe[0][0] + L2 * Xe[1][0] + L3 * Xe[2][0],
    //!     L1 * Xe[0][1] + L2 * Xe[1][1] + L3 * Xe[2][1]
    //! }; // ==> {0.2,0.3}
    //! @endcode
    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_tri3_2d(const MatXe& Xe, const VecP& p, Real& LI, Real& LJ, const int ix = 0, const int iy = 1)
    {
        // NOTE: we only use the x- y-coordinate
        //
        // Xe = {A,B,C}
        //
        static constexpr int I = 0;
        static constexpr int J = 1;
        static constexpr int K = 2;

        // area of I-J-p
        auto ax = Xe[J][ix] - Xe[I][ix];
        auto ay = Xe[J][iy] - Xe[I][iy];
        auto bx = p[ix] - Xe[J][ix];
        auto by = p[iy] - Xe[J][iy];
        auto sK = ax * by - ay * bx;
        // area of J-K-p
        ax = p[ix] - Xe[K][ix];
        ay = p[iy] - Xe[K][iy];
        auto sI = ay * bx - ax * by;
        // area of K-I-p
        bx = p[ix] - Xe[I][ix];
        by = p[iy] - Xe[I][iy];
        auto sJ = ax * by - ay * bx;

        auto ssum = sI + sJ + sK;
        LI = sI / ssum;
        LJ = sJ / ssum;
    }
    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_tri3_3d(const MatXe& Xe, const VecP& p, Real& LI, Real& LJ)
    {
        using namespace VecMatAlg;

        auto& xI = view3(&Xe[0][0]);
        auto& xJ = view3(&Xe[1][0]);
        auto& xK = view3(&Xe[2][0]);
        auto dIK = xI - xK;
        auto dJK = xJ - xK;
        auto d = xK - view3(&p[0]);
        auto a11 = dot(dIK, dIK);
        auto a12 = dot(dIK, dJK);
        auto a22 = dot(dJK, dJK);
        auto b1 = dot(dIK, d);
        auto b2 = dot(dJK, d);
        auto det = a11 * a22 - a12 * a12; _ASSERT(det != 0);
        LI = -( a22 * b1 - a12 * b2)/ det;
        LJ = -(-a12 * b1 + a11 * b2)/ det;
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_tri6_2d(const MatXe& Xe, const VecP& p, Real& LI, Real& LJ)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>; // std::remove_cvref_t<decltype(Xe[0][0])>
        
        static constexpr int NN = 6;
        static constexpr int ND = 2;
        static constexpr int NC = 2;

        type x[NC] = { type{1} / type{3},type{1} / type{3} };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_tri6(x[0], x[1], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_tri6(x[0], x[1], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        LI = x[0];
        LJ = x[1];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_tri6_3d(const MatXe& Xe, const VecP& p, Real& LI, Real& LJ)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 6;
        static constexpr int ND = 3;
        static constexpr int NC = 2;

        type x[NC] = { type{1} / type{3},type{1} / type{3} };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_tri6(x[0], x[1], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_tri6(x[0], x[1], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        LI = x[0];
        LJ = x[1];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_quad4_2d(const MatXe& Xe, const VecP& p, Real& s, Real& t)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 4;
        static constexpr int ND = 2;
        static constexpr int NC = 2;

        type x[NC] = { 0,0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_quad4(x[0], x[1], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_quad4(x[0], x[1], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        s = x[0];
        t = x[1];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_quad4_2d_analytic(const MatXe& Xe, const VecP& p, Real& s, Real& t, int ix = 0, int iy = 1)
    {
        //? This algorithm can get right results for degenerated element.

        static constexpr auto eps = std::numeric_limits<Real>::epsilon();

        // L         K
        // +---------+
        // |         |
        // |         |
        // |         |
        // +---------+
        // I         J

        auto xI = Xe[0][ix];      auto yI = Xe[0][iy];
        auto xJ = Xe[1][ix] - xI; auto yJ = Xe[1][iy] - yI;
        auto xK = Xe[2][ix] - xI; auto yK = Xe[2][iy] - yI;
        auto xL = Xe[3][ix] - xI; auto yL = Xe[3][iy] - yI;
        auto x  = p[ix] - xI;     auto y  = p[iy] - yI;

        auto a1 = xK - xJ - xL;
        auto b1 = yK - yJ - yL;
        auto d1 = xL * yJ - xJ * yL;

        // parallelogram
        if (std::fabs(a1) <= eps && std::fabs(b1) <= eps) {
            s = (2 * y * xL - 2 * x * yL) / d1 - Real{ 1 };
            t = (2 * x * yJ - 2 * y * xJ) / d1 - Real{ 1 };
        }
        else {
            auto d2 = xK * yJ - xJ * yK;
            auto ds = d1 - d2;
            // edge-IJ is parallel with edge-LK
            if (std::fabs(ds) <= eps) {
                t = (2 * x * yJ - 2 * y * xJ) / d1 - Real{ 1 };
                auto axJ = xK - xL;
                auto ayJ = yK - yL;
                auto det = y * (axJ - xJ) + x * (yJ - ayJ) - d1;
                s = (2 * x * yL - 2 * y * xL) / det - Real{ 1 };
                return;
            }
            auto d3 = xL * yK - xK * yL;
            auto dt = d1 - d3;
            // edge-IL is parallel with edge-JK
            if (std::fabs(dt) <= eps) {
                s = (2 * y * xL - 2 * x * yL) / d1 - Real{ 1 };
                auto bxL = xK - xJ;
                auto byL = yK - yJ;
                auto det = y * (xL - bxL) - x * (yL - byL) - d1;
                t = (2 * y * xJ - 2 * x * yJ) / det - Real{ 1 };
                return;
            }

            // general quadrilateral
            auto a = x * b1 - y * a1;
            auto e = a + d2;
            auto tmp = std::sqrt(e * e + ds * (
                2 * x * (-yJ + yK + yL) -
                2 * y * (-xJ + xK + xL) + d1 + d2));
            auto b = d2 + a;
            auto c = d3 - a;
            auto s1 = (b + tmp) / ds;
            auto s2 = (b - tmp) / ds;
            auto t1 = (tmp + c) / dt;
            auto t2 = (tmp - c) / dt;
            if (s1 < -1 || s1 > 1 || t1 < -1 || t1 > 1) {
                s = s2; t = t2;
            }
            else {
                s = s1; t = t1;
            }
        }
    }
        
    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_quad4_3d(const MatXe& Xe, const VecP& p, Real& s, Real& t)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 4;
        static constexpr int ND = 3;
        static constexpr int NC = 2;

        type x[NC] = { 0,0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_quad4(x[0], x[1], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_quad4(x[0], x[1], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        s = x[0];
        t = x[1];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_quad8_2d(const MatXe& Xe, const VecP& p, Real& s, Real& t)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 8;
        static constexpr int ND = 2;
        static constexpr int NC = 2;

        type x[NC] = { 0,0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_quad8(x[0], x[1], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_quad8(x[0], x[1], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        s = x[0];
        t = x[1];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_quad8_3d(const MatXe& Xe, const VecP& p, Real& s, Real& t)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;

        static constexpr int NN = 8;
        static constexpr int ND = 3;
        static constexpr int NC = 2;

        type x[NC] = { 0,0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_quad8(x[0], x[1], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_quad8(x[0], x[1], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        s = x[0];
        t = x[1];
    }

    //! @brief Inverse isoparametric mappings for tetrahedron
    //! @param Xe Vertex array for TET: {I,J,K,L}
    //! @param p  Query point: {x,y,z}
    //! @param L1 Local coordinate for I
    //! @param L2 Local coordinate for J
    //! @param L3 Local coordinate for K
    //! @code
    //! double Xe[4][3] = {
    //!     { 0,0,0 },
    //!     { 0,0,1 },
    //!     { 1,0,0 },
    //!     { 0,1,0 }};
    //! double p[3] = {0.1,0.2,0.3};
    //! double L1, L2, L3;
    //! inverse_iso_mapping_tet4(Xe, p, L1, L2, L3);
    //! // ==> L1=0.4, L2=0.3, L3=0.1
    //! double L4 = 1 - L1 - L2 - L3;
    //! double q[3] = {
    //!     L1* XeT[0][0] + L2 * XeT[1][0] + L3 * XeT[2][0] + L4 * XeT[3][0],
    //!     L1* XeT[0][1] + L2 * XeT[1][1] + L3 * XeT[2][1] + L4 * XeT[3][1],
    //!     L1* XeT[0][2] + L2 * XeT[1][2] + L3 * XeT[2][2] + L4 * XeT[3][2]
    //! };
    //! // ==> q = {0.1,0.2,0.3}
    //! @endcode
    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_tet4(const MatXe& Xe, const VecP& p, Real& L1, Real& L2, Real& L3)
    {
        static constexpr int I = 0;
        static constexpr int J = 1;
        static constexpr int K = 2;
        static constexpr int L = 3;

        // volume of Tetrahedron: 1/6*dot(a,cross(b,c))
        // volume of I-J-K-p
        auto ax = p[0] - Xe[I][0];
        auto ay = p[1] - Xe[I][1];
        auto az = p[2] - Xe[I][2];
        auto bx = Xe[J][0] - Xe[I][0]; // I-->J
        auto by = Xe[J][1] - Xe[I][1];
        auto bz = Xe[J][2] - Xe[I][2];
        auto cx = Xe[K][0] - Xe[J][0]; // J-->K
        auto cy = Xe[K][1] - Xe[J][1];
        auto cz = Xe[K][2] - Xe[J][2];
        auto sx = by * cz - bz * cy; // = cross(b,c)
        auto sy = bz * cx - bx * cz;
        auto sz = bx * cy - by * cx;
        auto vL = ax * sx + ay * sy + az * sz;

        // volume of I-L-J-p
        cx = Xe[I][0] - Xe[L][0]; // L-->I
        cy = Xe[I][1] - Xe[L][1];
        cz = Xe[I][2] - Xe[L][2];
        sx = by * cz - bz * cy; // = cross(b,c)
        sy = bz * cx - bx * cz;
        sz = bx * cy - by * cx;
        auto vK = ax * sx + ay * sy + az * sz;

        // volume of I-K-L-p
        bx = Xe[I][0] - Xe[K][0];
        by = Xe[I][1] - Xe[K][1];
        bz = Xe[I][2] - Xe[K][2];
        sx = by * cz - bz * cy; // = cross(b,c)
        sy = bz * cx - bx * cz;
        sz = bx * cy - by * cx;
        auto vJ = ax * sx + ay * sy + az * sz;

        // volume of J-L-K-p
        ax = p[0] - Xe[L][0];
        ay = p[1] - Xe[L][1];
        az = p[2] - Xe[L][2];
        bx = Xe[J][0] - Xe[K][0];
        by = Xe[J][1] - Xe[K][1];
        bz = Xe[J][2] - Xe[K][2];
        cx = Xe[L][0] - Xe[J][0];
        cy = Xe[L][1] - Xe[J][1];
        cz = Xe[L][2] - Xe[J][2];
        sx = by * cz - bz * cy; // = cross(b,c)
        sy = bz * cx - bx * cz;
        sz = bx * cy - by * cx;
        auto vI = ax * sx + ay * sy + az * sz;

        auto vSum = vI + vJ + vK + vL;
        L1 = vI / vSum;
        L2 = vJ / vSum;
        L3 = vK / vSum;
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_tet10(const MatXe& Xe, const VecP& p, Real& L1, Real& L2, Real& L3)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 10;
        static constexpr int ND = 3;
        static constexpr int NC = 3;

        type x[NC] = { Real{ 1 } / Real{ 4 } ,Real{ 1 } / Real{ 4 } , Real{ 1 } / Real{ 4 } };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_tet10(x[0], x[1], x[2], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_tet10(x[0], x[1], x[2], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        L1 = x[0];
        L2 = x[1];
        L3 = x[2];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_hex8(const MatXe& Xe, const VecP& p, Real& s, Real& t, Real& r)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 8;
        static constexpr int ND = 3;
        static constexpr int NC = 3;

        type x[NC] = { 0,0,0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_hex8(x[0], x[1], x[2], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_hex8(x[0], x[1], x[2], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        s = x[0];
        t = x[1];
        r = x[2];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_hex20(const MatXe& Xe, const VecP& p, Real& s, Real& t, Real& r)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>; // std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 20;
        static constexpr int ND = 3;
        static constexpr int NC = 3;

        type x[NC] = { 0,0,0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_hex20(x[0], x[1], x[2], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_hex20(x[0], x[1], x[2], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        s = x[0];
        t = x[1];
        r = x[2];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_pyramid5(const MatXe& Xe, const VecP& p, Real& s, Real& t, Real& r)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 5;
        static constexpr int ND = 3;
        static constexpr int NC = 3;

        type x[NC] = { 0,0,0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_pyramid5(x[0], x[1], x[2], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_pyramid5(x[0], x[1], x[2], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        s = x[0];
        t = x[1];
        r = x[2];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_pyramid13(const MatXe& Xe, const VecP& p, Real& s, Real& t, Real& r)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>
        
        static constexpr int NN = 13;
        static constexpr int ND = 3;
        static constexpr int NC = 3;

        type x[NC] = { 0,0,0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_pyramid13(x[0], x[1], x[2], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_pyramid13(x[0], x[1], x[2], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        s = x[0];
        t = x[1];
        r = x[2];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_prism6(const MatXe& Xe, const VecP& p, Real& L1, Real& L2, Real& r)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 6;
        static constexpr int ND = 3;
        static constexpr int NC = 3;

        type x[NC] = { type{ 1 } / type{ 3 } ,type{ 1 } / type{ 3 } , 0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_prism6(x[0], x[1], x[2], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_prism6(x[0], x[1], x[2], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        L1 = x[0];
        L2 = x[1];
        r  = x[2];
    }

    template<typename MatXe, typename VecP, typename Real>
    void inverse_iso_mapping_prism15(const MatXe& Xe, const VecP& p, Real& L1, Real& L2, Real& r)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(Xe[0][0])>>;// std::remove_cvref_t<decltype(Xe[0][0])>

        static constexpr int NN = 15;
        static constexpr int ND = 3;
        static constexpr int NC = 3;

        type x[NC] = { type{ 1 } / type{ 3 }, type{ 1 } / type{ 3 }, 0 };
        MinimizeLM<NC, ND>::solve(
            [&Xe](const tvec<type, NC>& x, tvec<type, ND>& f) {
                tvec<type, NN> fn;
                shape_function_prism15(x[0], x[1], x[2], fn);
                f = dot(fn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            [&Xe](const tvec<type, NC>& x, tmat<type, NC, ND>& J) {
                tmat<type, NC, NN> dfn;
                shape_function_diff_prism15(x[0], x[1], x[2], dfn);
                J = dot(dfn, mat_view<type, NN, ND>(&Xe[0][0]));
            },
            p,
            x
        );
        L1 = x[0];
        L2 = x[1];
        r = x[2];
    }
}
