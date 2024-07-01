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

#pragma once
//!-----------------------------------------------------------------------------
//! \file ShapeFunctions.hpp
//! \brief Define shape functions for standard 1-D/2-D/3-D finite elements.
//! \author ZHANG Bing, zhangbing@hfut.edu.cn
//! \date 2024-04-18
//! \copyright (c)2024 All rights reserved.
//!-----------------------------------------------------------------------------

namespace EasyLib {

    template<typename Real, typename Vec2>
    void shape_function_bar2(Real s, Vec2& N)
    {
        N[0] = Real{ 0.5 } *(Real{ 1 } - s);
        N[1] = Real{ 0.5 } *(Real{ 1 } + s);
    }
    template<typename Real, typename Vec2>
    void shape_function_diff_bar2(Real s, Vec2& dNds)
    {
        dNds[0] = -1;
        dNds[1] =  1;
    }

    template<typename Real, typename Vec3>
    void shape_function_bar3(Real s, Vec3& N)
    {
        N[0] = (s * s - s) * Real{ 0.5 };
        N[1] = (s * s + s) * Real{ 0.5 };
        N[2] = Real{ 1 } - s * s;
    }
    template<typename Real, typename Vec3>
    void shape_function_diff_bar3(Real s, Vec3& dNds)
    {
        dNds[0] = s - Real{ 0.5 };
        dNds[1] = s + Real{ 0.5 };
        dNds[2] = Real{ -2 } *s;
    }

    template<typename Real, typename Vec3>
    void shape_function_tri3(Real l1, Real l2, Vec3& N)
    {
        // Robert D. Concepts and applications of finite element analysis (4th). p260, 7.1-1
        //N[0] = Real{ 1 } - l1 - l2;
        //N[1] = l1;
        //N[2] = l2;

        N[0] = l1;
        N[1] = l2;
        N[2] = Real{ 1 } - l1 - l2;
    }
    template<typename Real, typename Mat2x3>
    void shape_function_diff_tri3(Real l1, Real l2, Mat2x3& dNdl)
    {
        // Robert D. Concepts and applications of finite element analysis (4th). p260, 7.1-1
        //dNdl[0][0] =-1; dNdl[0][1] = 1; dNdl[0][2] = 0;
        //dNdl[1][0] =-1; dNdl[1][1] = 0; dNdl[1][2] = 1;

        dNdl[0][0] = 1; dNdl[0][1] = 0; dNdl[0][2] = -1;
        dNdl[1][0] = 0; dNdl[1][1] = 1; dNdl[1][2] = -1;
    }

    template<typename Real, typename Vec6>
    void shape_function_tri6(Real l1, Real l2, Vec6& N)
    {
        const auto l3 = Real{ 1 } - l1 - l2;

        N[0] = l1 * (Real{ 2 } * l1 - Real{ 1 });
        N[1] = l2 * (Real{ 2 } * l2 - Real{ 1 });
        N[2] = l3 * (Real{ 2 } * l3 - Real{ 1 });

        N[3] = Real{ 4 } * l1 * l2;
        N[4] = Real{ 4 } * l2 * l3;
        N[5] = Real{ 4 } * l3 * l1;
    }
    template<typename Real, typename Mat2x6>
    void shape_function_diff_tri6(Real l1, Real l2, Mat2x6& dNdl)
    {
        const auto l3 = Real{ 1 } - l1 - l2;

        dNdl[0][0] =  4 * l1 - 1;
        dNdl[0][1] =  0;
        dNdl[0][2] = -4 * l3 + 1;
        dNdl[0][3] =  4 * l2;
        dNdl[0][4] = -4 * l2;
        dNdl[0][5] =  4 * (l3 - l1);

        dNdl[1][0] =  0;
        dNdl[1][1] =  4 * l2 - 1;
        dNdl[1][2] = -4 * l3 + 1;
        dNdl[1][3] =  4 * l1;
        dNdl[1][4] =  4 * (l3 - l2);
        dNdl[1][5] = -4 * l1;
    }

    template<typename Real, typename Vec4>
    void shape_function_quad4(Real s, Real t, Vec4& N)
    {
        N[0] = Real{ 0.25 } * (1 - s) * (1 - t);
        N[1] = Real{ 0.25 } * (1 + s) * (1 - t);
        N[2] = Real{ 0.25 } * (1 + s) * (1 + t);
        N[3] = Real{ 0.25 } * (1 - s) * (1 + t);
    }
    template<typename Real, typename Mat2x4>
    void shape_function_diff_quad4(Real s, Real t, Mat2x4& dN)
    {
        dN[0][0] = Real{ 0.25 } * ( t - 1);
        dN[0][1] = Real{ 0.25 } * ( 1 - t);
        dN[0][2] = Real{ 0.25 } * ( 1 + t);
        dN[0][3] = Real{ 0.25 } * (-1 - t);

        dN[1][0] = Real{ 0.25 } * ( s - 1);
        dN[1][1] = Real{ 0.25 } * (-1 - s);
        dN[1][2] = Real{ 0.25 } * ( 1 + s);
        dN[1][3] = Real{ 0.25 } * ( 1 - s);
    }

    template<typename Real, typename Vec8>
    void shape_function_quad8(Real s, Real t, Vec8& N)
    {
        N[0] = Real{ 0.25 } * (1 - s) * (1 - t) * (-s - t - 1);
        N[1] = Real{ 0.25 } * (1 + s) * (1 - t) * ( s - t - 1);
        N[2] = Real{ 0.25 } * (1 + s) * (1 + t) * ( s + t - 1);
        N[3] = Real{ 0.25 } * (1 - s) * (1 + t) * (-s + t - 1);
        N[4] = Real{ 0.50 } * (1 - s * s) * (1 - t);
        N[5] = Real{ 0.50 } * (1 - t * t) * (1 + s);
        N[6] = Real{ 0.50 } * (1 - s * s) * (1 + t);
        N[7] = Real{ 0.50 } * (1 - t * t) * (1 - s);
    }
    template<typename Real, typename Mat2x8>
    void shape_function_diff_quad8(Real s, Real t, Mat2x8& dN)
    {
        //0.25*(1-t)*(2*s+t),0.25*(1-t)*(2*s-t),0.25*(1+t)*(2*s+t),0.25*(1+t)*(2*s-t),s*(t-1),0.5*(1-t*t),-s*(t+1),0.5*(t*t-1)
        dN[0][0] = Real{ 0.25 } * (1 - t) * (2 * s + t);
        dN[0][1] = Real{ 0.25 } * (1 - t) * (2 * s - t);
        dN[0][2] = Real{ 0.25 } * (1 + t) * (2 * s + t);
        dN[0][3] = Real{ 0.25 } * (1 + t) * (2 * s - t);
        dN[0][4] =  s * (t - 1);
        dN[0][5] = Real{ 0.50 } * (1 - t * t);
        dN[0][6] = -s * (t + 1);
        dN[0][7] = Real{ 0.50 } * (t * t - 1);

        //0.25*(1-s)*(2*t+s),0.25*(1+s)*(2*t-s),0.25*(1+s)*(2*t+s),0.25*(1-s)*(2*t-s),0.5*(s*s-1),-t*(1+s),0.5*(1-s*s),t*(s-1)]
        dN[1][0] = Real{ 0.25 } * (1 - s) * (2 * t + s);
        dN[1][1] = Real{ 0.25 } * (1 + s) * (2 * t - s);
        dN[1][2] = Real{ 0.25 } * (1 + s) * (2 * t + s);
        dN[1][3] = Real{ 0.25 } * (1 - s) * (2 * t - s);
        dN[1][4] = Real{ 0.50 } * (s * s - 1);
        dN[1][5] = -t * (1 + s);
        dN[1][6] = Real{ 0.50 } * (1 - s * s);
        dN[1][7] =  t * (s - 1);
    }

    template<typename Real, typename Vec4>
    void shape_function_tet4(Real L1, Real L2, Real L3, Vec4& N)
    {
        N[0] = L1;
        N[1] = L2;
        N[2] = L3;
        N[3] = 1 - L1 - L2 - L3;
    }
    template<typename Real, typename Mat3x4>
    void shape_function_diff_tet4(Real L1, Real L2, Real L3, Mat3x4& dN)
    {
        dN[0][0] = 1; dN[0][1] = 0; dN[0][2] = 0; dN[0][3] = -1;
        dN[1][0] = 0; dN[1][1] = 1; dN[1][2] = 0; dN[1][3] = -1;
        dN[1][0] = 0; dN[1][1] = 0; dN[1][2] = 1; dN[1][3] = -1;
    }

    template<typename Real, typename Vec10>
    void shape_function_tet10(Real L1, Real L2, Real L3, Vec10& N)
    {
        const auto L4 = 1 - L1 - L2 - L3;

        N[0] = (2 * L1 - 1) * L1;
        N[1] = (2 * L2 - 1) * L2;
        N[2] = (2 * L3 - 1) * L3;
        N[3] = (2 * L4 - 1) * L4;
        N[4] = 4 * L1 * L2;
        N[5] = 4 * L2 * L3;
        N[6] = 4 * L1 * L3;
        N[7] = 4 * L1 * L4;
        N[8] = 4 * L2 * L4;
        N[9] = 4 * L3 * L4;
    }
    template<typename Real, typename Mat3x10>
    void shape_function_diff_tet10(Real L1, Real L2, Real L3, Mat3x10& dN)
    {
        const auto L4 = 1 - L1 - L2 - L3;

        dN[0][0] = -1 + 4 * L1;
        dN[0][1] = 0;
        dN[0][2] = 0;
        dN[0][3] = 1 - 4 * L4;
        dN[0][4] = 4 * L2;
        dN[0][5] = 0;
        dN[0][6] = 4 * L3;
        dN[0][7] = 4 * (L4 - L1);
        dN[0][8] = -4 * L2;
        dN[0][9] = -4 * L3;

        dN[1][0] = 0;
        dN[1][1] = -1 + 4 * L2;
        dN[1][2] = 0;
        dN[1][3] = 1 - 4 * L4;
        dN[1][4] = 4 * L1;
        dN[1][5] = 4 * L3;
        dN[1][6] = 0;
        dN[1][7] = -4 * L1;
        dN[1][8] = 4 * (L4 - L2);
        dN[1][9] = -4 * L3;

        dN[2][0] = 0;
        dN[2][1] = 0;
        dN[2][2] = -1 + 4 * L3;
        dN[2][3] = 1 - 4 * L4;
        dN[2][4] = 0;
        dN[2][5] = 4 * L2;
        dN[2][6] = 4 * L1;
        dN[2][7] = -4 * L1;
        dN[2][8] = -4 * L2;
        dN[2][9] = 4 * (L4 - L3);
    }

    template<typename Real, typename Vec5>
    void shape_function_pyramid5(Real s, Real t, Real r, Vec5& N)
    {
        N[0] = Real{ 0.125 } * (1 - s) * (1 - t) * (1 - r);
        N[1] = Real{ 0.125 } * (1 + s) * (1 - t) * (1 - r);
        N[2] = Real{ 0.125 } * (1 + s) * (1 + t) * (1 - r);
        N[3] = Real{ 0.125 } * (1 - s) * (1 + t) * (1 - r);
        N[4] = Real{ 0.500 } * (1 + r);
    }
    template<typename Real, typename Mat3x5>
    void shape_function_diff_pyramid5(Real s, Real t, Real r, Mat3x5& dN)
    {
        dN[0][0] = Real{ -0.125 } *(-1 + r) * (-1 + t);
        dN[0][1] = Real{  0.125 } *(-1 + r) * (-1 + t);
        dN[0][2] = Real{ -0.125 } *(-1 + r) * ( 1 + t);
        dN[0][3] = Real{  0.125 } *(-1 + r) * ( 1 + t);
        dN[0][4] = 0;

        dN[1][0] = Real{ -0.125 } *(-1 + r) * (-1 + s);
        dN[1][1] = Real{  0.125 } *(-1 + r) * ( 1 + s);
        dN[1][2] = Real{ -0.125 } *(-1 + r) * ( 1 + s);
        dN[1][3] = Real{  0.125 } *(-1 + r) * (-1 + s);
        dN[1][4] = 0;

        dN[2][0] = Real{ -0.125 } *(-1 + s) * (-1 + t);
        dN[2][1] = Real{  0.125 } *(1 + s) * (-1 + t);
        dN[2][2] = Real{ -0.125 } *(1 + s) * (1 + t);
        dN[2][3] = Real{  0.125 } *(-1 + s) * (1 + t);
        dN[2][4] = Real{  0.500 };
    }

    template<typename Real, typename Vec13>
    void shape_function_pyramid13(Real s, Real t, Real r, Vec13& N)
    {
        const auto q = Real{ 0.5 } * (1 - r);
        N[0] = Real{ 0.25 } * q * (1 - s) * (1 - t) * (-1 - q * s - q * t);
        N[1] = Real{ 0.25 } * q * (1 + s) * (1 - t) * (-1 + q * s - q * t);
        N[2] = Real{ 0.25 } * q * (1 + s) * (1 + t) * (-1 + q * s + q * t);
        N[3] = Real{ 0.25 } * q * (1 - s) * (1 + t) * (-1 - q * s + q * t);
        N[4] = (1 - q) * (1 - 2 * q);

        N[5] = Real{ 0.50 } * q * q * (1 - t) * (1 - s * s);
        N[6] = Real{ 0.50 } * q * q * (1 + s) * (1 - t * t);
        N[7] = Real{ 0.50 } * q * q * (1 + t) * (1 - s * s);
        N[8] = Real{ 0.50 } * q * q * (1 - s) * (1 - t * t);

        const auto p = q * (1 - q);
        N[ 9] = p * (1 - s - t + s * t);
        N[10] = p * (1 + s - t - s * t);
        N[11] = p * (1 + s + t + s * t);
        N[12] = p * (1 - s + t - s * t);
    }
    template<typename Real, typename Mat3x13>
    void shape_function_diff_pyramid13(Real s, Real t, Real r, Mat3x13& dN)
    {
        auto rm1 = r - 1;
        auto rp1 = r + 1;
        auto rm1_rm1 = rm1 * rm1;
        auto rm1_rp1 = rm1 * rp1;
        auto rm1_tp1 = rm1 * ( 1 + t);
        auto rm1_tm1 = rm1 * (-1 + t);
        auto rm1_sm1 = rm1 * (-1 + s);
        auto rm1_sp1 = rm1 * ( 1 + s);
        auto sm1_tm1 = (-1 + s) * (-1 + t);
        auto sm1_tp1 = (-1 + s) * ( 1 + t);
        auto sp1_tm1 = ( 1 + s) * (-1 + t);
        auto sp1_tp1 = ( 1 + s) * ( 1 + t);
        auto rm1s = rm1 * s;
        auto rm1t = rm1 * t;
        dN[0][ 0] = Real{ 0.0625 } * rm1_tm1 * ( rp1 - 2 * rm1s - rm1t);
        dN[0][ 1] = Real{ 0.0625 } * rm1_tm1 * (-rp1 - 2 * rm1s + rm1t);
        dN[0][ 2] = Real{ 0.0625 } * rm1_tp1 * ( rp1 + 2 * rm1s + rm1t);
        dN[0][ 3] = Real{ 0.0625 } * rm1_tp1 * (-rp1 + 2 * rm1s - rm1t);
        dN[0][ 4] = 0;
        dN[0][ 5] = Real{  0.250 } *rm1_rm1 * (-1 + t) * s;
        dN[0][ 6] = Real{ -0.125 } *rm1_rm1 * (-1 + t * t);
        dN[0][ 7] = Real{ -0.250 } *rm1_rm1 * ( 1 + t) * s;
        dN[0][ 8] = Real{  0.125 } *rm1_rm1 * (-1 + t * t);
        dN[0][ 9] = Real{ -0.250 } *rm1_rp1 * (-1 + t);
        dN[0][10] = Real{  0.250 } *rm1_rp1 * (-1 + t);
        dN[0][11] = Real{ -0.250 } *rm1_rp1 * ( 1 + t);
        dN[0][12] = Real{  0.250 } *rm1_rp1 * ( 1 + t);

        dN[1][ 0] = Real{  0.0625 } *rm1_sm1 * ( rp1 - rm1s - 2 * rm1t);
        dN[1][ 1] = Real{  0.0625 } *rm1_sp1 * (-rp1 - rm1s + 2 * rm1t);
        dN[1][ 2] = Real{  0.0625 } *rm1_sp1 * ( rp1 + rm1s + 2 * rm1t);
        dN[1][ 3] = Real{  0.0625 } *rm1_sm1 * (-rp1 + rm1s - 2 * rm1t);
        dN[1][ 4] = 0;
        dN[1][ 5] = Real{  0.125 } * rm1_rm1 * (-1 + s * s);
        dN[1][ 6] = Real{ -0.250 } * rm1_rm1 * ( 1 + s) * t;
        dN[1][ 7] = Real{ -0.125 } * rm1_rm1 * (-1 + s * s);
        dN[1][ 8] = Real{  0.250 } * rm1_rm1 * (-1 + s) * t;
        dN[1][ 9] = Real{ -0.250 } * rm1_rp1 * (-1 + s);
        dN[1][10] = Real{  0.250 } * rm1_rp1 * ( 1 + s);
        dN[1][11] = Real{ -0.250 } * rm1_rp1 * ( 1 + s);
        dN[1][12] = Real{  0.250 } * rm1_rp1 * (-1 + s);

        dN[2][ 0] = Real{ -0.125 } *sm1_tm1 * (-1 + rm1 * (s + t));
        dN[2][ 1] = Real{ -0.125 } *sp1_tm1 * ( 1 + rm1 * (s - t));
        dN[2][ 2] = Real{  0.125 } *sp1_tp1 * ( 1 + rm1 * (s + t));
        dN[2][ 3] = Real{  0.125 } *sm1_tp1 * (-1 + rm1 * (s - t));
        dN[2][ 4] = Real{  0.500 } + r;
        dN[2][ 5] = Real{  0.250 } * rm1_tm1 * (-1 + s * s);
        dN[2][ 6] = Real{ -0.250 } * rm1_sp1 * (-1 + t * t);
        dN[2][ 7] = Real{ -0.250 } * rm1_tp1 * (-1 + s * s);
        dN[2][ 8] = Real{  0.250 } * rm1_sm1 * (-1 + t * t);
        dN[2][ 9] = Real{ -0.500 } * r * sm1_tm1;
        dN[2][10] = Real{  0.500 } * r * sp1_tm1;
        dN[2][11] = Real{ -0.500 } * r * sp1_tp1;
        dN[2][12] = Real{  0.500 } * r * sm1_tp1;
    }

    template<typename Real, typename Vec6>
    void shape_function_prism6(Real L1, Real L2, Real r, Vec6& N)
    {
        const auto L3 = 1 - L1 - L2;
        N[0] = Real{ 0.5 } *L1 * (1 - r);
        N[1] = Real{ 0.5 } *L2 * (1 - r);
        N[2] = Real{ 0.5 } *L3 * (1 - r);
        N[3] = Real{ 0.5 } *L1 * (1 + r);
        N[4] = Real{ 0.5 } *L2 * (1 + r);
        N[5] = Real{ 0.5 } *L3 * (1 + r);
    }
    template<typename Real, typename Mat3x6>
    void shape_function_diff_prism6(Real L1, Real L2, Real r, Mat3x6& dN)
    {
        dN[0][0] = Real{ 0.5 } * (1 - r);
        dN[0][1] = 0;
        dN[0][2] = Real{ 0.5 } * (-1 + r);
        dN[0][3] = Real{ 0.5 } * ( 1 + r);
        dN[0][4] = 0;
        dN[0][5] = Real{ 0.5 } * (-1 - r);

        dN[1][0] = 0;
        dN[1][1] = Real{ 0.5 } * ( 1 - r);
        dN[1][2] = Real{ 0.5 } * (-1 + r);
        dN[1][3] = 0;
        dN[1][4] = Real{ 0.5 } * ( 1 + r);
        dN[1][5] = Real{ 0.5 } * (-1 - r);

        const auto L3 = 1 - L1 - L2;
        dN[2][0] = Real{ -0.5 } * L1;
        dN[2][1] = Real{ -0.5 } * L2;
        dN[2][2] = Real{ -0.5 } * L3;
        dN[2][3] = Real{  0.5 } * L1;
        dN[2][4] = Real{  0.5 } * L2;
        dN[2][5] = Real{  0.5 } * L3;
    }

    template<typename Real, typename Vec15>
    void shape_function_prism15(Real L1, Real L2, Real r, Vec15& N)
    {
        const auto L3 = 1 - L1 - L2;

        N[0] = Real{ 0.5 } * L1 * (2 * L1 - 1) * (1 - r) - L1 * (1 - r * r);
        N[1] = Real{ 0.5 } * L2 * (2 * L2 - 1) * (1 - r) - L2 * (1 - r * r);
        N[2] = Real{ 0.5 } * L3 * (2 * L3 - 1) * (1 - r) - L3 * (1 - r * r);
        N[3] = Real{ 0.5 } * L1 * (2 * L1 - 1) * (1 + r) - L1 * (1 - r * r);
        N[4] = Real{ 0.5 } * L2 * (2 * L2 - 1) * (1 + r) - L2 * (1 - r * r);
        N[5] = Real{ 0.5 } * L3 * (2 * L3 - 1) * (1 + r) - L3 * (1 - r * r);

        N[6] = 2 * L1 * L2 * (1 - r);
        N[7] = 2 * L2 * L3 * (1 - r);
        N[8] = 2 * L3 * L1 * (1 - r);

        N[ 9] = 2 * L1 * L2 * (1 + r);
        N[10] = 2 * L2 * L3 * (1 + r);
        N[11] = 2 * L3 * L1 * (1 + r);

        N[12] = L1 * (1 - r * r);
        N[13] = L2 * (1 - r * r);
        N[14] = L3 * (1 - r * r);
    }
    template<typename Real, typename Mat3x15>
    void shape_function_diff_prism15(Real L1, Real L2, Real r, Mat3x15& dN)
    {
        const auto L3 = 1 - L1 - L2;

        dN[0][ 0] = Real{ -0.5 } * (-3 + 4 * L1 - 2 * r) * (-1 + r);
        dN[0][ 1] = 0;
        dN[0][ 2] = Real{ -0.5 } * ( 3 - 4 * L3 + 2 * r) * (-1 + r);
        dN[0][ 3] = Real{  0.5 } * (-3 + 4 * L1 + 2 * r) * (+1 + r);
        dN[0][ 4] = 0;
        dN[0][ 5] = Real{  0.5 } * ( 3 - 4 * L3 - 2 * r) * (+1 + r);
        dN[0][ 6] = 2 * L2 * (+1 - r);
        dN[0][ 7] = 2 * L2 * (-1 + r);
        dN[0][ 8] = 2 * (L1 - L3) * (-1 + r);
        dN[0][ 9] = 2 * L2 * (+1 + r);
        dN[0][10] = 2 * L2 * (-1 - r);
        dN[0][11] = 2 * (L1 - L3) * (-1 - r);
        dN[0][12] = +1 - r * r;
        dN[0][13] = 0;
        dN[0][14] = -1 + r * r;

        dN[1][ 0] = 0;
        dN[1][ 1] = Real{ -0.5 } * (-3 + 4 * L2 - 2 * r) * (-1 + r);
        dN[1][ 2] = Real{ -0.5 } * ( 3 - 4 * L3 + 2 * r) * (-1 + r);
        dN[1][ 3] = 0;
        dN[1][ 4] = Real{  0.5 } * (-3 + 4 * L2 + 2 * r) * (1 + r);
        dN[1][ 5] = Real{  0.5 } * ( 3 - 4 * L3 - 2 * r) * (1 + r);
        dN[1][ 6] = 2 * L1 * (+1 - r);
        dN[1][ 7] = 2 * (L2 - L3) * (-1 + r);
        dN[1][ 8] = 2 * L1 * (-1 + r);
        dN[1][ 9] = 2 * L1 * (+1 + r);
        dN[1][10] = 2 * (L2 - L3) * (-1 - r);
        dN[1][11] = 2 * L1 * (-1 - r);
        dN[1][12] = 0;
        dN[1][13] = +1 - r * r;
        dN[1][14] = -1 + r * r;

        dN[2][ 0] = Real{ 0.5 } *L1 * ( 1 - 2 * L1 + 4 * r);
        dN[2][ 1] = Real{ 0.5 } *L2 * ( 1 - 2 * L2 + 4 * r);
        dN[2][ 2] = Real{ 0.5 } *L3 * ( 1 - 2 * L3 + 4 * r);
        dN[2][ 3] = Real{ 0.5 } *L1 * (-1 + 2 * L1 + 4 * r);
        dN[2][ 4] = Real{ 0.5 } *L2 * (-1 + 2 * L2 + 4 * r);
        dN[2][ 5] = Real{ 0.5 } *L3 * (-1 + 2 * L3 + 4 * r);
        dN[2][ 6] = -2 * L1 * L2;
        dN[2][ 7] = -2 * L2 * L3;
        dN[2][ 8] = -2 * L1 * L3;
        dN[2][ 9] =  2 * L1 * L2;
        dN[2][10] =  2 * L2 * L3;
        dN[2][11] =  2 * L1 * L3;
        dN[2][12] = -2 * L1 * r;
        dN[2][13] = -2 * L2 * r;
        dN[2][14] = -2 * L3 * r;
    }

    template<typename Real, typename Vec8>
    void shape_function_hex8(Real s, Real t, Real r, Vec8& N)
    {
        N[0] = Real{ 0.125 } * (1 - s) * (1 - t) * (1 - r);
        N[1] = Real{ 0.125 } * (1 + s) * (1 - t) * (1 - r);
        N[2] = Real{ 0.125 } * (1 + s) * (1 + t) * (1 - r);
        N[3] = Real{ 0.125 } * (1 - s) * (1 + t) * (1 - r);
        N[4] = Real{ 0.125 } * (1 - s) * (1 - t) * (1 + r);
        N[5] = Real{ 0.125 } * (1 + s) * (1 - t) * (1 + r);
        N[6] = Real{ 0.125 } * (1 + s) * (1 + t) * (1 + r);
        N[7] = Real{ 0.125 } * (1 - s) * (1 + t) * (1 + r);
    }
    template<typename Real, typename Mat3x8>
    void shape_function_diff_hex8(Real s, Real t, Real r, Mat3x8& dN)
    {
        dN[0][0] = Real{ -0.125 } * (1 - t) * (1 - r);
        dN[0][1] = Real{  0.125 } * (1 - t) * (1 - r);
        dN[0][2] = Real{  0.125 } * (1 + t) * (1 - r);
        dN[0][3] = Real{ -0.125 } * (1 + t) * (1 - r);
        dN[0][4] = Real{ -0.125 } * (1 - t) * (1 + r);
        dN[0][5] = Real{  0.125 } * (1 - t) * (1 + r);
        dN[0][6] = Real{  0.125 } * (1 + t) * (1 + r);
        dN[0][7] = Real{ -0.125 } * (1 + t) * (1 + r);

        dN[1][0] = Real{ -0.125 } * (1 - s) * (1 - r);
        dN[1][1] = Real{ -0.125 } * (1 + s) * (1 - r);
        dN[1][2] = Real{  0.125 } * (1 + s) * (1 - r);
        dN[1][3] = Real{  0.125 } * (1 - s) * (1 - r);
        dN[1][4] = Real{ -0.125 } * (1 - s) * (1 + r);
        dN[1][5] = Real{ -0.125 } * (1 + s) * (1 + r);
        dN[1][6] = Real{  0.125 } * (1 + s) * (1 + r);
        dN[1][7] = Real{  0.125 } * (1 - s) * (1 + r);

        dN[2][0] = Real{ -0.125 } * (1 - s) * (1 - t);
        dN[2][1] = Real{ -0.125 } * (1 + s) * (1 - t);
        dN[2][2] = Real{ -0.125 } * (1 + s) * (1 + t);
        dN[2][3] = Real{ -0.125 } * (1 - s) * (1 + t);
        dN[2][4] = Real{  0.125 } * (1 - s) * (1 - t);
        dN[2][5] = Real{  0.125 } * (1 + s) * (1 - t);
        dN[2][6] = Real{  0.125 } * (1 + s) * (1 + t);
        dN[2][7] = Real{  0.125 } * (1 - s) * (1 + t);
    }

    template<typename Real, typename Vec20>
    void shape_function_hex20(Real s, Real t, Real r, Vec20& N)
    {
        N[0] = Real{ 0.125 } * (1 - s) * (1 - t) * (1 - r) * (-s - t - r - 2);
        N[1] = Real{ 0.125 } * (1 + s) * (1 - t) * (1 - r) * ( s - t - r - 2);
        N[2] = Real{ 0.125 } * (1 + s) * (1 + t) * (1 - r) * ( s + t - r - 2);
        N[3] = Real{ 0.125 } * (1 - s) * (1 + t) * (1 - r) * (-s + t - r - 2);

        N[4] = Real{ 0.125 } * (1 - s) * (1 - t) * (1 + r) * (-s - t + r - 2);
        N[5] = Real{ 0.125 } * (1 + s) * (1 - t) * (1 + r) * ( s - t + r - 2);
        N[6] = Real{ 0.125 } * (1 + s) * (1 + t) * (1 + r) * ( s + t + r - 2);
        N[7] = Real{ 0.125 } * (1 - s) * (1 + t) * (1 + r) * (-s + t + r - 2);

        N[ 8] = Real{ 0.25 } * (1 - s * s) * (1 - t) * (1 - r);
        N[ 9] = Real{ 0.25 } * (1 + s) * (1 - t * t) * (1 - r);
        N[10] = Real{ 0.25 } * (1 - s * s) * (1 + t) * (1 - r);
        N[11] = Real{ 0.25 } * (1 - s) * (1 - t * t) * (1 - r);

        N[12] = Real{ 0.25 } * (1 - s * s) * (1 - t) * (1 + r);
        N[13] = Real{ 0.25 } * (1 + s) * (1 - t * t) * (1 + r);
        N[14] = Real{ 0.25 } * (1 - s * s) * (1 + t) * (1 + r);
        N[15] = Real{ 0.25 } * (1 - s) * (1 - t * t) * (1 + r);

        N[16] = Real{ 0.25 } * (1 - s) * (1 - t) * (1 - r * r);
        N[17] = Real{ 0.25 } * (1 + s) * (1 - t) * (1 - r * r);
        N[18] = Real{ 0.25 } * (1 + s) * (1 + t) * (1 - r * r);
        N[19] = Real{ 0.25 } * (1 - s) * (1 + t) * (1 - r * r);
    }
    template<typename Real, typename Mat3x20>
    void shape_function_diff_hex20(Real s, Real t, Real r, Mat3x20& dN)
    {
        dN[0][ 0] = Real{  0.125 } * (1 - t) * (1 - r) * (1 + 2 * s + t + r);
        dN[0][ 1] = Real{ -0.125 } * (1 - t) * (1 - r) * (1 - 2 * s + t + r);
        dN[0][ 2] = Real{ -0.125 } * (1 + t) * (1 - r) * (1 - 2 * s - t + r);
        dN[0][ 3] = Real{  0.125 } * (1 + t) * (1 - r) * (1 + 2 * s - t + r);
        dN[0][ 4] = Real{  0.125 } * (1 - t) * (1 + r) * (1 + 2 * s + t - r);
        dN[0][ 5] = Real{ -0.125 } * (1 - t) * (1 + r) * (1 - 2 * s + t - r);
        dN[0][ 6] = Real{ -0.125 } * (1 + t) * (1 + r) * (1 - 2 * s - t - r);
        dN[0][ 7] = Real{  0.125 } * (1 + t) * (1 + r) * (1 + 2 * s - t - r);
        dN[0][ 8] = Real{ -0.500 } * s * (1 - t) * (1 - r);
        dN[0][ 9] = Real{  0.250 } * (1 - t * t) * (1 - r);
        dN[0][10] = Real{ -0.500 } * s * (1 + t) * (1 - r);
        dN[0][11] = Real{ -0.250 } * (1 - t * t) * (1 - r);
        dN[0][12] = Real{ -0.500 } * s * (1 - t) * (1 + r);
        dN[0][13] = Real{  0.250 } * (1 - t * t) * (1 + r);
        dN[0][14] = Real{ -0.500 } * s * (1 + t) * (1 + r);
        dN[0][15] = Real{ -0.250 } * (1 - t * t) * (1 + r);
        dN[0][16] = Real{ -0.250 } * (1 - t) * (1 - r * r);
        dN[0][17] = Real{  0.250 } * (1 - t) * (1 - r * r);
        dN[0][18] = Real{  0.250 } * (1 + t) * (1 - r * r);
        dN[0][19] = Real{ -0.250 } * (1 + t) * (1 - r * r);

        dN[1][ 0] = Real{  0.125 } * (1 - s) * (1 - r) * (1 + s + 2 * t + r);
        dN[1][ 1] = Real{  0.125 } * (1 + s) * (1 - r) * (1 - s + 2 * t + r);
        dN[1][ 2] = Real{ -0.125 } * (1 + s) * (1 - r) * (1 - s - 2 * t + r);
        dN[1][ 3] = Real{ -0.125 } * (1 - s) * (1 - r) * (1 + s - 2 * t + r);
        dN[1][ 4] = Real{  0.125 } * (1 - s) * (1 + r) * (1 + s + 2 * t - r);
        dN[1][ 5] = Real{  0.125 } * (1 + s) * (1 + r) * (1 - s + 2 * t - r);
        dN[1][ 6] = Real{ -0.125 } * (1 + s) * (1 + r) * (1 - s - 2 * t - r);
        dN[1][ 7] = Real{ -0.125 } * (1 - s) * (1 + r) * (1 + s - 2 * t - r);
        dN[1][ 8] = Real{ -0.250 } * (1 - s * s) * (1 - r);
        dN[1][ 9] = Real{ -0.500 } * (1 + s) * t * (1 - r);
        dN[1][10] = Real{  0.250 } * (1 - s * s) * (1 - r);
        dN[1][11] = Real{ -0.500 } * (1 - s) * t * (1 - r);
        dN[1][12] = Real{ -0.250 } * (1 - s * s) * (1 + r);
        dN[1][13] = Real{ -0.500 } * (1 + s) * t * (1 + r);
        dN[1][14] = Real{  0.250 } * (1 - s * s) * (1 + r);
        dN[1][15] = Real{ -0.500 } * (1 - s) * t * (1 + r);
        dN[1][16] = Real{ -0.250 } * (1 - s) * (1 - r * r);
        dN[1][17] = Real{ -0.250 } * (1 + s) * (1 - r * r);
        dN[1][18] = Real{  0.250 } * (1 + s) * (1 - r * r);
        dN[1][19] = Real{  0.250 } * (1 - s) * (1 - r * r);

        dN[2][ 0] = Real{  0.125 } * (1 - s) * (1 - t) * (1 + s + t + 2 * r);
        dN[2][ 1] = Real{  0.125 } * (1 + s) * (1 - t) * (1 - s + t + 2 * r);
        dN[2][ 2] = Real{  0.125 } * (1 + s) * (1 + t) * (1 - s - t + 2 * r);
        dN[2][ 3] = Real{  0.125 } * (1 - s) * (1 + t) * (1 + s - t + 2 * r);
        dN[2][ 4] = Real{ -0.125 } * (1 - s) * (1 - t) * (1 + s + t - 2 * r);
        dN[2][ 5] = Real{ -0.125 } * (1 + s) * (1 - t) * (1 - s + t - 2 * r);
        dN[2][ 6] = Real{ -0.125 } * (1 + s) * (1 + t) * (1 - s - t - 2 * r);
        dN[2][ 7] = Real{ -0.125 } * (1 - s) * (1 + t) * (1 + s - t - 2 * r);
        dN[2][ 8] = Real{ -0.250 } * (1 - s * s) * (1 - t);
        dN[2][ 9] = Real{ -0.250 } * (1 + s) * (1 - t * t);
        dN[2][10] = Real{ -0.250 } * (1 - s * s) * (1 + t);
        dN[2][11] = Real{ -0.250 } * (1 - s) * (1 - t * t);
        dN[2][12] = Real{  0.250 } * (1 - s * s) * (1 - t);
        dN[2][13] = Real{  0.250 } * (1 + s) * (1 - t * t);
        dN[2][14] = Real{  0.250 } * (1 - s * s) * (1 + t);
        dN[2][15] = Real{  0.250 } * (1 - s) * (1 - t * t);
        dN[2][16] = Real{ -0.500 } * (1 - s) * (1 - t) * r;
        dN[2][17] = Real{ -0.500 } * (1 + s) * (1 - t) * r;
        dN[2][18] = Real{ -0.500 } * (1 + s) * (1 + t) * r;
        dN[2][19] = Real{ -0.500 } * (1 - s) * (1 + t) * r;
    }
}
