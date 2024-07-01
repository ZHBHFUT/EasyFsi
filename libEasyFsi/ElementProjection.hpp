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
//! \file ElementProjection.hpp
//! \brief Do element projection for standard 2-D/3-D line and surface elements.
//! \author ZHANG Bing, zhangbing@hfut.edu.cn
//! \date 2024-04-18
//! \copyright (c)2024 All rights reserved.
//!-----------------------------------------------------------------------------
#pragma once
#include "VecMatAlg.hpp"
#include "ShapeFunctions.hpp"
#include "InverseIsoMapping.hpp"

namespace EasyLib {

    //-----------------------------------------------------------------------------
    //            compute element projection for standard elements
    //-----------------------------------------------------------------------------
    // Point in element:
    //  X(s,t) = {N1(s,t),N2(s,t),...}.{{xI,yI},{xJ,yJ},...}
    //         =          N           .       Xe
    // 
    //  for given p={x,y}, find {s,t}:
    //       min(dot(X({s,t}) - p), dot(X({s,t}) - p))
    // 
    //  def: F({s,t}) = dot(X({s,t}) - p, X({s,t}) - p)
    // 
    //  solve:
    //       dF/d{s,t} = 0
    // 
    // ==>  f({s,t}) = dot(dN . Xe, X({s,t}) - p) = 0
    // 
    // ==>
    //      J = df/d{s,t} = [J11,J12]
    //                      [J21,J22]
    // 
    //   J11 =       dot(ddN/dds . Xe, X({s,t}) - p) + dot(dNds . Xe, dNds . Xe)
    //   J12 = J21 = dot(ddN/dsdt. Xe, X({s,t}) - p) + dot(dNds . Xe, dNdt . Xe)
    //   J22 =       dot(ddN/ddt . Xe, X({s,t}) - p) + dot(dNdt . Xe, dNdt . Xe)
    //
    //   {s,t}_{k+1} = {s,t}_k - J^-1 . f({s,t}_k)
    //
    //-----------------------------------------------------------------------------
 
    //! projection point to line segment.
    //! @param [in]  Xe vertex of line segment, one row per point.
    //! @param [in]  p  point need to projection
    //! @param [out] s  output parameter, in [-1,1]
    //! @return true=succeed, false=nearest point is not on the line segment.
    //! @note The projection point = (1-s)/2*a + (1+s)/2*b
    template<typename Real, typename Mat2x2, typename Vec2>
    bool projection_bar2_2d(const Mat2x2& Xe, const Vec2& p, Real& s)
    {
        using namespace VecMatAlg;
        auto ab = view2(&Xe[1][0]) - view2(&Xe[0][0]);
        auto ap = view2(&p[0]) - view2(&Xe[0][0]);
        auto l2 = norm_sq(ab);
        s = l2 != 0 ? 2 * dot(ab, ap) / l2 - 1 : -2;
        if      (s < -1) { s = -1; return false; }
        else if (s > +1) { s = +1; return false; }
        else return true;
    }
    //! projection point to line segment.
    //! @param [in]  a  start vertex of line segment
    //! @param [in]  b  end vertex of line segment
    //! @param [in]  p  point need to projection
    //! @param [out] s  output parameter, in [-1,1]
    //! @return true=succeed, false=nearest point is not on the line segment.
    //! @note The projection point = (1-s)/2*a + (1+s)/2*b
    template<typename Real, typename Vec2>
    bool projection_bar2_2d(const Vec2& a, const Vec2& b, const Vec2& p, Real& s)
    {
        using namespace VecMatAlg;
        auto ab = view2(&b[0]) - view2(&a[0]);
        auto ap = view2(&p[0]) - view2(&a[0]);
        auto l2 = norm_sq(ab);
        s = l2 != 0 ? 2 * dot(ab, ap) / l2 - 1 : -2;
        if      (s < -1) { s = -1; return false; }
        else if (s > +1) { s = +1; return false; }
        else return true;
    }

    template<typename Real, typename Mat2x3, typename Vec3>
    bool projection_bar2_3d(const Mat2x3& Xe, const Vec3& p, Real& s)
    {
        using namespace VecMatAlg;

        auto ab = view3(&Xe[1][0]) - view3(&Xe[0][0]);
        auto ap = view3(&p[0]) - view3(&Xe[0][0]);
        auto l2 = norm_sq(ab);
        s = l2 != 0 ? 2 * dot(ab, ap) / l2 - 1 : -2;
        if      (s < -1) { s = -1; return false; }
        else if (s > +1) { s = +1; return false; }
        else return true;
    }
    template<typename Real, typename Vec3>
    bool projection_bar2_3d(const Vec3& a, const Vec3& b, const Vec3& p, Real& s)
    {
        using namespace VecMatAlg;

        auto ab = view3(&b[0]) - view3(&a[0]);
        auto ap = view3(&p[0]) - view3(&a[0]);
        auto l2 = norm_sq(ab);
        s = l2 != 0 ? 2 * dot(ab, ap) / l2 - 1 : -2;
        if      (s < -1) { s = -1; return false; }
        else if (s > +1) { s = +1; return false; }
        else return true;
    }

    //! @brief projection point to line segment.
    //! @param Xe 3*2 matrix, vertex of line segment, one row per point.
    //! @param p  point need to projection
    //! @param s  output parameter, in [-1,1]
    //! @return true=succeed, false=nearest point is not on the line segment.
    //! @code
    //! double Xe[3][2] = { 0,0,1,0,0.5,0.2 };
    //! double p[2] = { 0.0,0.5 };
    //! double s;
    //! projection_bar3_2d(Xe, p, s); // s=-0.62331444575555510
    //! assert(std::fabs(s + 0.62331444575555510) <= std::numeric_limits<double>::epsilon())
    //! double fn[3];
    //! shape_function_bar3(s, fn);
    //! auto q = dot(mat_view<double, 1, NN>(fn), mat_view<double, NN, ND>(&Xe[0][0]));
    //! auto d = distance(vec_view<double, ND>(p), q);
    //! assert(std::fabs(d - 0.42205858482545345) <= std::numeric_limits<double>::epsilon())
    //! @endcode
    template<typename Real, typename Mat3x2, typename Vec2>
    bool projection_bar3_2d(const Mat3x2& Xe, const Vec2& p, Real& s)
    {
        inverse_iso_mapping_bar3_2d(Xe, p, s);

        // fix
        if      (s < -1.0) { s = -1.0; return false; }
        else if (s > +1.0) { s = +1.0; return false; }
        return true;
    }
    template<typename Real, typename Vec2>
    bool projection_bar3_2d(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& p, Real& s)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(a[0])>>;

        tvec<type, 2> Xe[3] = {
            {a[0],a[1]},
            {b[0],b[1]},
            {c[0],c[1]}
        };

        inverse_iso_mapping_bar3_2d(Xe, p, s);

        // fix
        if      (s < -1.0) { s = -1.0; return false; }
        else if (s > +1.0) { s = +1.0; return false; }
        return true;
    }

    //! @brief projection point to line segment.
    //! @param Xe 3*3 matrix, vertex of line segment, one row per point.
    //! @param p  size=3, point need to projection
    //! @param s  output parameter, in [-1,1]
    //! @return true=succeed, false=nearest point is not on the line segment.
    //! @code
    //! double Xe[3][3] = { 0,0,0,1,0,0,0.5,0.2,0.2 };
    //! double p[3] = { 0.1,0.0,0.0 };
    //! double s;
    //! projection_bar3_3d(Xe, p, s); // s = -0.90490722258298975
    //! assert(std::fabs(s + 0.90490722258298975) <= std::numeric_limits<double>::epsilon());
    //! double fn[3];
    //! shape_function_bar3(s, fn);
    //! auto q = dot(mat_view<double, 1, 3>(fn), mat_view<double, 3, 2>(&Xe[0][0]));
    //! auto d = distance(vec_view<double, 3>(p), q);
    //! _ASSERT(std::fabs(d - 0.073323951692689070) <= std::numeric_limits<double>::epsilon());
    //! @endcode
    template<typename Real, typename Mat3x3, typename Vec3>
    bool projection_bar3_3d(const Mat3x3& Xe, const Vec3& p, Real& s)
    {
        inverse_iso_mapping_bar3_3d(Xe, p, s);

        // fix
        if      (s < -1.0) { s = -1.0; return false; }
        else if (s > +1.0) { s = +1.0; return false; }
        return true;
    }
    template<typename Real, typename Vec3>
    bool projection_bar3_3d(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& p, Real& s)
    {
        using namespace VecMatAlg;
        using type = std::remove_cv_t<std::remove_reference_t<decltype(a[0])>>;

        tvec<type, 3> Xe[3] = {
            {a[0],a[1],a[2]},
            {b[0],b[1],b[2]},
            {c[0],c[1],c[2]}
        };

        inverse_iso_mapping_bar3_3d(Xe, p, s);

        // fix
        if      (s < -1.0) { s = -1.0; return false; }
        else if (s >  1.0) { s =  1.0; return false; }
        return true;
    }

    //! projection point to triangle.
    //! @param [in]  Xe   3*3 vertex matrix, one row for a vertex
    //! @param [in]  p    point need to project
    //! @param [out] l1   area coordinate, in [0,1]
    //! @param [out] l2   area coordinate, in [0,1]
    //! @return true=succeed, false=nearest point is not in the triangle.
    //! @note The projection point = l1*Xe(0,:) + l2*Xe(1,:) + (1-l1-l2)*Xe(2,:)
    //! @code
    //! constexpr int NN = 3;
    //! constexpr int ND = 3;
    //! double Xe[NN][ND] = { 0,0,0,1,0,0,0.5,0.5,0 };
    //! double p[ND] = { 0,0.5,0 };
    //! double l1, l2;
    //! projection_tri3(Xe, p, l1, l2); // => l1=0.5, l2=0
    //! double fn[NN];
    //! shape_function_tri3(l1, l2, fn);
    //! auto q = dot(mat_view<double, 1, NN>(fn), mat_view<double, NN, ND>(&Xe[0][0]));
    //! auto d = distance(vec_view<double, ND>(p), q); // => d = sqrt(2)/4
    //! _ASSERT(std::fabs(l1 - 0.5) < std::numeric_limits<double>::epsilon());
    //! _ASSERT(std::fabs(l2) < std::numeric_limits<double>::epsilon());
    //! _ASSERT(std::fabs(d - 0.25 * sqrt(2)) < std::numeric_limits<double>::epsilon());
    //! @endcode
    template<typename Real, typename Mat3x3, typename Vec3>
    bool projection_tri3(const Mat3x3& Xe, const Vec3& p, Real& l1, Real& l2)
    {
        using namespace VecMatAlg;
        
        inverse_iso_mapping_tri3_3d(Xe, p, l1, l2);

        auto l3 = 1 - l1 - l2;

        // triangle
        if (l1 >= 0 && l2 >= 0 && l3 >= 0) {
            return true;
        }
        else if (l3 <= 0) {
            if (l2 <= 0) {
                l1 = 1; l2 = 0;
            }
            else if (l1 <= 0) {
                l1 = 0; l2 = 1;
            }
            else {
                // project to line(v0,v1)
                projection_bar2_3d(view3(&Xe[0][0]), view3(&Xe[1][0]), view3(&p[0]), l2);//l2:[-1,1]
                l2 = Real{ 0.5 } * (l2 + 1.0);//to [0,1]
                l1 = 1 - l2;
            }
        }
        else if (l1 <= 0) {
            if (l2 <= 0) {
                l1 = l2 = 0;
            }
            else {
                // project to line(v1,v2)
                projection_bar2_3d(view3(&Xe[1][0]), view3(&Xe[2][0]), view3(&p[0]), l3);//l3:[-1,1]
                l3 = Real{ 0.5 } * (l3 + 1.0);//to [0,1]
                l2 = 1 - l3;
                l1 = 0;
            }
        }
        else {
            // project to line(v0,v2)
            projection_bar2_3d(view3(&Xe[0][0]), view3(&Xe[2][0]), view3(&p[0]), l3);//l3:[-1,1]
            l3 = Real{ 0.5 } * (l3 + 1.0);//to [0,1]
            l1 = 1 - l3;
            l2 = 0;
        }
        return false;
    }

    //! projection point to triangle.
    //! @param [in]  Xe   6*3 vertex matrix, one row for a vertex
    //! @param [in]  p    point need to project
    //! @param [out] l1   area coordinate, in [0,1]
    //! @param [out] l2   area coordinate, in [0,1]
    //! @return true=succeed, false=nearest point is not in the triangle.
    //! @code
    //! constexpr int NN = 6;
    //! constexpr int ND = 3;
    //! double Xe[NN][ND] = { 0,0,0, 1,0,0, 0.5,1,0, 0.5,0,0.25, 0.75,0.5,0.25, 0.25,0.5,0.25};
    //! double p[ND] = { 0.5,0.5,1.0 };
    //! double l1, l2;
    //! projection_tri6(Xe, p, l1, l2); // l1=l2=0.29182754901359119
    //! double fn[NN];
    //! shape_function_tri6(l1, l2, fn);
    //! auto q = dot(mat_view<double, 1, NN>(fn), mat_view<double, NN, ND>(&Xe[0][0]));
    //! auto d = distance(vec_view<double, ND>(p), q); // d=0.67702307980624926
    //! _ASSERT(std::fabs(l1 - 0.29182754901359119) < std::numeric_limits<double>::epsilon());
    //! _ASSERT(std::fabs(l2 - 0.29182754901359119) < std::numeric_limits<double>::epsilon());
    //! _ASSERT(std::fabs(d - 0.67702307980624926) < std::numeric_limits<double>::epsilon());
    //! @endcode
    template<typename Real, typename Mat6x3, typename Vec3>
    bool projection_tri6(const Mat6x3& Xe, const Vec3& p, Real& l1, Real& l2)
    {
        using namespace VecMatAlg;
        
        inverse_iso_mapping_tri6_3d(Xe, p, l1, l2);

        auto l3 = 1 - l1 - l2;

        //TRI
        if (l1 >= 0 && l2 >= 0 && l3 >= 0) {
            //do nothing
            return true;
        }
        else if (l3 <= 0) {
            if (l2 <= 0) {
                l1 = 1; l2 = 0;
            }
            else if (l1 <= 0) {
                l1 = 0; l2 = 1;
            }
            else {
                //project to line(v0,v1,v3)
                projection_bar3_3d(view3(&Xe[0][0]), view3(&Xe[1][0]), view3(&Xe[3][0]), view3(&p[0]), l2);//l2:[-1,1]
                l2 = Real{ 0.5 } * (l2 + 1);//to [0,1]
                l1 = 1 - l2;
            }
        }
        else if (l1 <= 0) {
            if (l2 <= 0) {
                l1 = l2 = 0;
            }
            else {
                //project to line(v1,v2,v4)
                projection_bar3_3d(view3(&Xe[1][0]), view3(&Xe[2][0]), view3(&Xe[4][0]), view3(&p[0]), l3);//l3:[-1,1]
                l3 = Real{ 0.5 } * (l3 + 1);//to [0,1]
                l2 = 1 - l3;
                l1 = 0;
            }
        }
        else {
            //project to line(v0,v2,v5)
            projection_bar3_3d(view3(&Xe[0][0]), view3(&Xe[2][0]), view3(&Xe[5][0]), view3(&p[0]), l3);//l3:[-1,1]
            l3 = Real{ 0.5 } * (l3 + 1);//to [0,1]
            l1 = 1 - l3;
            l2 = 0;
        }
        return false;
    }

    //! @brief projection point to quadrilateral.
    //! @param Xe  4*3 vertex matrix, one row for a vertex
    //! @param p   point need to project
    //! @param s   local coordinate, in [-1,1]
    //! @param t   local coordinate, in [-1,1]
    //! @return true=succeed, false=nearest point is not in the triangle.
    //! @code
    //! constexpr int NN = 4;
    //! constexpr int ND = 3;
    //! double Xe[NN][ND] = { 0,0,0, 1,0,0, 1,1,0.1, 0,1,0 };
    //! double p[ND] = { 0.25,0.25,0.25 };
    //! double s, t;
    //! projection_quad4(Xe, p, s, t); // s=t=-0.48752459474145472
    //! double fn[NN];
    //! shape_function_quad4(s, t, fn);
    //! auto q = dot(mat_view<double, 1, NN>(fn), mat_view<double, NN, ND>(&Xe[0][0]));
    //! auto d = distance(vec_view<double, ND>(p), q); // d=0.24359400499715687
    //! _ASSERT(std::fabs(s + 0.48752459474145472) < std::numeric_limits<double>::epsilon());
    //! _ASSERT(std::fabs(t + 0.48752459474145472) < std::numeric_limits<double>::epsilon());
    //! _ASSERT(std::fabs(d - 0.24359400499715687) < std::numeric_limits<double>::epsilon());
    //! @endcode
    template<typename Real, typename Mat4x3, typename Vec3>
    bool projection_quad4(const Mat4x3& Xe, const Vec3& p, Real& s, Real& t)
    {
        using namespace VecMatAlg;
        
        inverse_iso_mapping_quad4_3d(Xe, p, s, t);

        if      (s < -1) {
            projection_bar2_3d(view3(&Xe[0][0]), view3(&Xe[3][0]), view3(&p[0]), t);
            s = -1;
            return false;
        }
        else if (s > +1) {
            projection_bar2_3d(view3(&Xe[1][0]), view3(&Xe[2][0]), view3(&p[0]), t);
            s = 1;
            return false;
        }
        else if (t < -1) {
            projection_bar2_3d(view3(&Xe[0][0]), view3(&Xe[1][0]), view3(&p[0]), s);
            t = -1;
            return false;
        }
        else if (t > +1) {
            projection_bar2_3d(view3(&Xe[3][0]), view3(&Xe[2][0]), view3(&p[0]), s);
            t = 1;
            return false;
        }
        return true;
    }

    //! @brief projection point to quadrilateral.
    //! @param Xe  8*3 vertex matrix, one row for a vertex
    //! @param p   point need to project
    //! @param s   local coordinate, in [-1,1]
    //! @param t   local coordinate, in [-1,1]
    //! @return true=succeed, false=nearest point is not in the triangle.
    //! @code
    //! constexpr int NN = 8;
    //! constexpr int ND = 3;
    //! double Xe[NN][ND] = { 0,0,0, 1,0,0, 1,1,0, 0,1,0, 0.5,0,0.2, 1,0.5,0.2, 0.5,1,0.2, 0,0.5,0.2 };
    //! double p[ND] = { 0.25,0.25,0.5 };
    //! double s, t;
    //! projection_quad8(Xe, p, s, t); // s=t=-0.39661342471961217
    //! double fn[NN];
    //! shape_function_quad8(s, t, fn);
    //! auto q = dot(mat_view<double, 1, NN>(fn), mat_view<double, NN, ND>(&Xe[0][0]));
    //! auto d = distance(vec_view<double, ND>(p), q); // d=0.17857101176790449
    //! _ASSERT(std::fabs(s + 0.39661342471961217) < std::numeric_limits<double>::epsilon());
    //! _ASSERT(std::fabs(t + 0.39661342471961217) < std::numeric_limits<double>::epsilon());
    //! _ASSERT(std::fabs(d - 0.17857101176790449) < std::numeric_limits<double>::epsilon());
    //! @endcode
    template<typename Real, typename Mat8x3, typename Vec3>
    bool projection_quad8(const Mat8x3& Xe, const Vec3& p, Real& s, Real& t)
    {
        using namespace VecMatAlg;
        
        inverse_iso_mapping_quad8_3d(Xe, p, s, t);

        if      (s < -1) {
            projection_bar3_3d(view3(&Xe[0][0]), view3(&Xe[3][0]), view3(&Xe[7][0]), view3(&p[0]), t);
            s = -1;
            return false;
        }
        else if (s > +1) {
            projection_bar3_3d(view3(&Xe[1][0]), view3(&Xe[2][0]), view3(&Xe[5][0]), view3(&p[0]), t);
            s = 1;
            return false;
        }
        else if (t < -1) {
            projection_bar3_3d(view3(&Xe[0][0]), view3(&Xe[1][0]), view3(&Xe[4][0]), view3(&p[0]), s);
            t = -1;
            return false;
        }
        else if (t > +1) {
            projection_bar3_3d(view3(&Xe[3][0]), view3(&Xe[2][0]), view3(&Xe[6][0]), view3(&p[0]), s);
            t = 1;
            return false;
        }
        return true;
    }

}
