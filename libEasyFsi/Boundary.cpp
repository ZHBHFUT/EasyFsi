#include <cmath>
#include <numbers>
#include <span>

#include "Logger.hpp"
#include "LinAlgs.h"
#include "ProjectToElement.hpp"
#include "Boundary.hpp"

namespace EasyLib {
    const int npf[FT_POLYGON + 1] = {
        2, // FT_BAR2 = 0,
        3, // FT_BAR3,
        3, // FT_TRI3,
        6, // FT_TRI6,
        4, // FT_QUAD4,
        8, // FT_QUAD8,
        0  // FT_POLYGON
    };
    
    const int face_order[FT_POLYGON + 1] = {
        1, // FT_BAR2 = 0,
        2, // FT_BAR3,
        1, // FT_TRI3,
        2, // FT_TRI6,
        1, // FT_QUAD4,
        2, // FT_QUAD8,
        1  // FT_POLYGON
    };
    static const ZoneTopo ft2zt[FT_POLYGON + 1] = {
        ZT_CURVE,   // FT_BAR2 = 0,
        ZT_CURVE,   // FT_BAR3,
        ZT_SURFACE, // FT_TRI3,
        ZT_SURFACE, // FT_TRI6,
        ZT_SURFACE, // FT_QUAD4,
        ZT_SURFACE, // FT_QUAD8,
        ZT_SURFACE  // FT_POLYGON
    };

    static void compute_fc(int nf, const Boundary::vec3* pnts, const FaceTopo* ftopo, const MeshConnectivity& face_nodes, Boundary::vec3* fcent)
    {
        double fn[npf_max];
        for (int i = 0; i < nf; ++i) {
            auto nodes = face_nodes[i];
            if      (ftopo[i] == FT_BAR2) {
                auto& x0 = pnts[nodes[0]];
                auto& x1 = pnts[nodes[1]];
                fcent[i] = 0.5 * (x0 + x1);
            }
            else if (ftopo[i] == FT_BAR3) {
                fcent[i] = pnts[nodes[2]];
            }
            else if (ftopo[i] == FT_TRI3) {
                auto& x0 = pnts[nodes[0]];
                auto& x1 = pnts[nodes[1]];
                auto& x2 = pnts[nodes[2]];
                fcent[i] = (x0 + x1 + x2) / 3.0;
            }
            else if (ftopo[i] == FT_TRI6) {
                fn[0] = -1.0 / 9.0;
                fn[1] = -1.0 / 9.0;
                fn[2] = -1.0 / 9.0;
                fn[3] = 4.0 / 9.0;
                fn[4] = 4.0 / 9.0;
                fn[5] = 4.0 / 9.0;
                auto& x0 = pnts[nodes[0]];
                auto& x1 = pnts[nodes[1]];
                auto& x2 = pnts[nodes[2]];
                auto& x3 = pnts[nodes[3]];
                auto& x4 = pnts[nodes[4]];
                auto& x5 = pnts[nodes[5]];
                fcent[i].x = fn[0] * x0.x + fn[1] * x1.x + fn[2] * x2.x + fn[3] * x3.x + fn[4] * x4.x + fn[5] * x5.x;
                fcent[i].y = fn[0] * x0.y + fn[1] * x1.y + fn[2] * x2.y + fn[3] * x3.y + fn[4] * x4.y + fn[5] * x5.y;
                fcent[i].z = fn[0] * x0.z + fn[1] * x1.z + fn[2] * x2.z + fn[3] * x3.z + fn[4] * x4.z + fn[5] * x5.z;
            }
            else if (ftopo[i] == FT_QUAD4) {
                //  3-----2
                //  |   / |
                //  | /   |
                //  0-----1
                auto& x0 = pnts[nodes[0]];
                auto& x1 = pnts[nodes[1]];
                auto& x2 = pnts[nodes[2]];
                auto& x3 = pnts[nodes[3]];
                auto s1 = cross(x1 - x0, x2 - x1).norm(); // 0-1-2
                auto s2 = cross(x2 - x0, x3 - x2).norm(); // 0-2-3
                auto s3 = cross(x1 - x0, x3 - x1).norm(); // 0-1-3
                auto s4 = cross(x2 - x1, x3 - x2).norm(); // 1-2-3
                auto c1 = (x0 + x1 + x2) / 3.0;
                auto c2 = (x0 + x2 + x3) / 3.0;
                auto c3 = (x0 + x1 + x3) / 3.0;
                auto c4 = (x1 + x2 + x3) / 3.0;
                auto ss = s1 + s2 + s3 + s4;
                fcent[i] = s1 / ss * c1;
                fcent[i] += s2 / ss * c2;
                fcent[i] += s3 / ss * c3;
                fcent[i] += s4 / ss * c4;
            }
            else if (ftopo[i] == FT_QUAD8) {
                fn[0] = -0.25;
                fn[1] = -0.25;
                fn[2] = -0.25;
                fn[3] = -0.25;
                fn[4] = 0.5;
                fn[5] = 0.5;
                fn[6] = 0.5;
                fn[7] = 0.5;
                auto& x0 = pnts[nodes[0]];
                auto& x1 = pnts[nodes[1]];
                auto& x2 = pnts[nodes[2]];
                auto& x3 = pnts[nodes[3]];
                auto& x4 = pnts[nodes[4]];
                auto& x5 = pnts[nodes[5]];
                auto& x6 = pnts[nodes[6]];
                auto& x7 = pnts[nodes[7]];
                fcent[i].x = fn[0] * x0.x + fn[1] * x1.x + fn[2] * x2.x + fn[3] * x3.x + fn[4] * x4.x + fn[5] * x5.x + fn[6] * x6.x + fn[7] * x7.x;
                fcent[i].y = fn[0] * x0.y + fn[1] * x1.y + fn[2] * x2.y + fn[3] * x3.y + fn[4] * x4.y + fn[5] * x5.y + fn[6] * x6.y + fn[7] * x7.y;
                fcent[i].z = fn[0] * x0.z + fn[1] * x1.z + fn[2] * x2.z + fn[3] * x3.z + fn[4] * x4.z + fn[5] * x5.z + fn[6] * x6.z + fn[7] * x7.z;
            }
            else if (ftopo[i] == FT_POLYGON) {
                auto& x0 = pnts[nodes[0]];
                double ss = 0;
                auto& fc = fcent[i];
                fc.assign(0, 0, 0);
                for (size_t j = 1; j < nodes.size() - 1; ++j) {
                    auto& x1 = pnts[nodes[j]];
                    auto& x2 = pnts[nodes[j + 1]];
                    auto s = cross(x1 - x0, x2 - x1).norm();
                    auto c = (x0 + x1 + x2) / 3.0;
                    ss += s;
                    fc += s * c;
                }
                fc /= ss;
            }
        }
    }

    static void compute_farea(int nf, const Boundary::vec3* pnts, const FaceTopo* ftopo, const MeshConnectivity& face_nodes, Boundary::vec3* farea)
    {
        //double fn[8];
        for (int i = 0; i < nf; ++i) {
            auto nodes = face_nodes[i];
            if (ftopo[i] == FT_BAR2 || ftopo[i] == FT_BAR3) {
                auto& x0 = pnts[nodes[0]];
                auto& x1 = pnts[nodes[1]];
                farea[i] = cross(x1 - x0, Boundary::vec3{ 0,0,1 });
            }
            //else if (ftopo[i] == FT_BAR3) {
            //    ShapeFunctionLINE3(0.0, fn);
            //    auto& x0 = pnts[nodes[0]];
            //    auto& x1 = pnts[nodes[1]];
            //    auto& x2 = pnts[nodes[2]];
            //    fcent[i] = fn[0] * x0 + fn[1] * x1 + fn[2] * x2;
            //}
            else if (ftopo[i] == FT_TRI3 || ftopo[i] == FT_TRI6) {
                auto& x0 = pnts[nodes[0]];
                auto& x1 = pnts[nodes[1]];
                auto& x2 = pnts[nodes[2]];
                farea[i] = cross(x1 - x0, x2 - x1) / 2.0;
            }
            //else if (ftopo[i] == FT_TRI6) {
            //    ShapeFunctionTRI6(1.0 / 3.0, 1.0 / 3.0, fn);
            //    auto& x0 = pnts[nodes[0]];
            //    auto& x1 = pnts[nodes[1]];
            //    auto& x2 = pnts[nodes[2]];
            //    auto& x3 = pnts[nodes[3]];
            //    auto& x4 = pnts[nodes[4]];
            //    auto& x5 = pnts[nodes[5]];
            //    fcent[i].x = fn[0] * x0.x + fn[1] * x1.x + fn[2] * x2.x + fn[3] * x3.x + fn[4] * x4.x + fn[5] * x5.x;
            //    fcent[i].y = fn[0] * x0.y + fn[1] * x1.y + fn[2] * x2.y + fn[3] * x3.y + fn[4] * x4.y + fn[5] * x5.y;
            //    fcent[i].z = fn[0] * x0.z + fn[1] * x1.z + fn[2] * x2.z + fn[3] * x3.z + fn[4] * x4.z + fn[5] * x5.z;
            //}
            else if (ftopo[i] == FT_QUAD4 || ftopo[i] == FT_QUAD8) {
                //  3-----2
                //  |   / |
                //  | /   |
                //  0-----1
                auto& x0 = pnts[nodes[0]];
                auto& x1 = pnts[nodes[1]];
                auto& x2 = pnts[nodes[2]];
                auto& x3 = pnts[nodes[3]];
                auto s1 = 0.5 * cross(x1 - x0, x2 - x1); // 0-1-2
                auto s2 = 0.5 * cross(x2 - x0, x3 - x2); // 0-2-3
                auto s3 = 0.5 * cross(x1 - x0, x3 - x1); // 0-1-3
                auto s4 = 0.5 * cross(x2 - x1, x3 - x2); // 1-2-3
                auto ss = s1 + s2 + s3 + s4;
                farea[i] = ss * 0.5;
            }
            //else if (ftopo[i] == FT_QUAD8) {
            //    ShapeFunctionQUAD8(0.0, 0.0, fn);
            //    auto& x0 = pnts[nodes[0]];
            //    auto& x1 = pnts[nodes[1]];
            //    auto& x2 = pnts[nodes[2]];
            //    auto& x3 = pnts[nodes[3]];
            //    auto& x4 = pnts[nodes[4]];
            //    auto& x5 = pnts[nodes[5]];
            //    auto& x6 = pnts[nodes[6]];
            //    auto& x7 = pnts[nodes[7]];
            //    fcent[i].x = fn[0] * x0.x + fn[1] * x1.x + fn[2] * x2.x + fn[3] * x3.x + fn[4] * x4.x + fn[5] * x5.x + fn[6] * x6.x + fn[7] * x7.x;
            //    fcent[i].y = fn[0] * x0.y + fn[1] * x1.y + fn[2] * x2.y + fn[3] * x3.y + fn[4] * x4.y + fn[5] * x5.y + fn[6] * x6.y + fn[7] * x7.y;
            //    fcent[i].z = fn[0] * x0.z + fn[1] * x1.z + fn[2] * x2.z + fn[3] * x3.z + fn[4] * x4.z + fn[5] * x5.z + fn[6] * x6.z + fn[7] * x7.z;
            //}
            else if (ftopo[i] == FT_POLYGON) {
                auto& x0 = pnts[nodes[0]];
                auto& fs = farea[i];
                fs.assign(0, 0, 0);
                for (size_t j = 1; j < nodes.size() - 1; ++j) {
                    auto& x1 = pnts[nodes[j + 0]];
                    auto& x2 = pnts[nodes[j + 1]];
                    auto s = 0.5 * cross(x1 - x0, x2 - x1);
                    fs += s;
                }
            }
        }
    }

    //! @brief 计算给定点集合的几何拓扑形式和局部坐标系
    //! @param np        点个数
    //! @param xyz       点坐标数组，存储格式为：[x0,y0,z0, x1,y1,z0, ..., xn,yn,zn]
    //! @param mt4x4_g2l 将原坐标系点变换到局部坐标系的变换矩阵，size=4*4，按行优先存储.
    //! @param biased_angle_deg 判断共面的最大偏斜角度，0表示严格共面.
    //! @return 返回拓扑结构类型，0=点(所有点重合)，1=共线(colinear)，2=共面(coplaner)，3=一般三维分布
    //! @note 当拓扑形式为一般三维分布时输出的变换矩阵为单位矩阵.
    static ZoneShape compute_topo_3d(int np, const double* xyz, double* mt4x4_g2l, double biased_angle_deg = 5)
    {
        // 算法流程：
        // 1. 计算局部坐标系原点 origin：选取第一个点
        // 2. 查找与原点最远的点作为 x 轴上的点 px，需要遍历所有点
        // 3. 根据 origin、px 计算 x 轴 axis_x：
        // 4. 查找与原点、px构成三角形面积最大的点 q，需要遍历所有点
        // 5. 根据origin、px、q计算局部坐标系 z 轴 axis_z
        // 6. 计算局部坐标系y轴：axis_y
        // 7. 形成坐标变换矩阵

        //
        // 变换矩阵格式：
        //          ( ax.x  ax.y  ax.z  -dot(ax,o) )
        //  T_g2l = ( ay.x  ay.y  ay.z  -dot(ay,o) )
        //          ( az.x  az.y  az.z  -dot(az,o) )
        //          (    0     0     0     1       )
        // 其中，o为局部坐标系原点坐标，ax、ay、az为局部坐标系三个坐标轴在原坐标系中的单位方向矢量
        // 
        // 于是，原坐标系点可通过如下公式变换到局部坐标系： Local = T_g2l * Global
        //

        using vec3 = Boundary::vec3;
        using real = double;

        static constexpr auto eps = std::numeric_limits<real>::epsilon();
        static constexpr auto rmax = std::numeric_limits<real>::max();
        const auto cos_err = std::abs(std::cos(std::numbers::pi_v<real> *(90 - biased_angle_deg) / 180));

        //--- 0. 初始化矩阵为单位矩阵
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                mt4x4_g2l[i * 4 + j] = i != j ? 0 : 1;

        if (np <= 1) { return ZS_POINT; }

        std::span<const vec3> pnts(reinterpret_cast<const vec3*>(xyz), np);

        auto& axis_x = vec3::view(mt4x4_g2l, 0); // x轴
        auto& axis_y = vec3::view(mt4x4_g2l, 4); // y轴
        auto& axis_z = vec3::view(mt4x4_g2l, 8); // z轴

        //--- 1. 局部坐标系原点：取第一个点

        auto& origin = pnts.front();

        //--- 2. 计算x轴：查找与原点最远的点作为局部坐标系x轴上的点

        vec3 px = pnts[1];
        {
            real dmax = distance_sq(px, origin);
            for (auto& c : pnts) {
                auto d = (distance_sq(c, origin));
                if (d > dmax) {
                    px = c;
                    dmax = d;
                }
            }
            if (dmax < 1E-30)return ZS_POINT;
        }
        axis_x = px - origin;
        // 单位化x轴，如果模值为0则表明所有点重合
        if (axis_x.normalize() <= eps) {
            axis_x.assign(1, 0, 0);
            return ZS_POINT;
        }

        //--- 3. 计算与原点、x轴最远点构成的三角形面积最大的点

        vec3 q = px;
        {
            real smax = 0;
            for (auto& c : pnts) {
                auto s = cross(axis_x, c - px).norm_sq();
                if (s > smax) {
                    q = c;
                    smax = s;
                }
            }
        }

        //--- 4. 计算局部坐标系的z轴：垂直于 o-px-q
        axis_z = cross(axis_x, q - origin);
        // 单位化z轴，如果模值为 0 则表明所有点共线
        if (axis_z.normalize() <= eps) {
            // 选取一个与x轴不相同的矢量作为z轴
            axis_z = cross(axis_x, vec3{ 1,0,0 });
            auto szmax = axis_z.norm_sq();
            auto tz = cross(axis_x, vec3{ 0,1,0 });
            auto sz = tz.norm_sq();
            if (sz > szmax) { szmax = sz; axis_z = tz; }
            tz = cross(axis_x, vec3{ 0,1,0 });
            sz = tz.norm_sq();
            if (sz > szmax) { szmax = sz; axis_z = tz; }

            // 计算y轴
            axis_y = cross(axis_z, axis_x);

            //--- 5. 更新局部坐标系变换矩阵的第4列

            mt4x4_g2l[3] = -dot(axis_x, origin);
            mt4x4_g2l[7] = -dot(axis_y, origin);
            mt4x4_g2l[11] = -dot(axis_z, origin);

            return ZS_COLINEAR;
        }
        // 其它情况，检查是否共面
        else {
            // 计算y轴
            axis_y = cross(axis_z, axis_x);

            // 计算所有点与 x-y 平面的夹角
            for (auto& c : pnts) {
                auto p = c - origin;
                p.normalize();

                // 夹角超过最大容差，表明点非共面
                if (std::abs(dot(axis_z, p)) > cos_err) {
                    axis_x.assign(1, 0, 0);
                    axis_y.assign(0, 1, 0);
                    axis_z.assign(0, 0, 1);
                    return ZS_GENERAL;
                }
            }

            // 一些特殊情形

            auto zx = std::abs(dot(axis_z, vec3{ 1,0,0 }));
            auto zy = std::abs(dot(axis_z, vec3{ 0,1,0 }));
            auto zz = std::abs(dot(axis_z, vec3{ 0,0,1 }));

            // x-y plane
            if (zx < 1E-8 && zy < 1E-8) {
                axis_x.assign(1, 0, 0);
                axis_y.assign(0, 1, 0);
                axis_z.assign(0, 0, 1);
            }
            // y-z plane
            else if (zy < 1E-8 && zz < 1E-8) {
                axis_x.assign(0, 1, 0);
                axis_y.assign(0, 0, 1);
                axis_z.assign(1, 0, 0);
            }
            // x-z plane
            else if (zx < 1E-8 && zz < 1E-8) {
                axis_x.assign(1, 0, 0);
                axis_y.assign(0, 0, 1);
                axis_z.assign(0, -1, 0);
            }

            //--- 5. 更新局部坐标系变换矩阵的第4列

            mt4x4_g2l[3] = -dot(axis_x, origin);
            mt4x4_g2l[7] = -dot(axis_y, origin);
            mt4x4_g2l[11] = -dot(axis_z, origin);

            return ZS_COPLANER;
        }
    }

    //! @brief 计算给定平面点集合的几何拓扑形式和局部坐标系
    //! @param np        点个数
    //! @param xy        点坐标数组，存储格式为：[x0,y0, x1,y1, ..., xn,yn]
    //! @param mt4x4_g2l 将原坐标系点变换到局部坐标系的变换矩阵，size=3*3，按行优先存储.
    //! @param biased_angle_deg 判断共面的最大偏斜角度，0表示严格共面.
    //! @return 返回拓扑结构类型，0=点(所有点重合)，1=共线(colinear)，2=共面(coplaner)
    //! @note 当拓扑形式为一般二维分布时输出的变换矩阵为单位矩阵.
    //static int compute_topo_2d(int np, const double* xy, double* mt3x3_g2l, double biased_angle_deg = 5)
    //{
    //    using vec2 = TinyVector<double, 2>;
    //
    //    //
    //    // a . b = |a| * |b| * cos(theta)
    //    //
    //    //  ==> cos(theta) = (a . b)/(|a|*|b|)
    //
    //    const auto cos_err = std::abs(std::cos(std::numbers::pi_v<double> *(biased_angle_deg) / 180));
    //
    //    mt3x3_g2l[0] = 1; mt3x3_g2l[1] = 0; mt3x3_g2l[2] = 0;
    //    mt3x3_g2l[3] = 0; mt3x3_g2l[4] = 1; mt3x3_g2l[5] = 0;
    //    mt3x3_g2l[6] = 0; mt3x3_g2l[7] = 0; mt3x3_g2l[8] = 1;
    //
    //    //--- 1. 查找与第1个点的最远点，作为向量 a
    //
    //    int i, n1 = 1;
    //    double dmax = 0;
    //    for (i = 1; i < np; ++i) {
    //        const auto dx = xy[2 * i + 0] - xy[0];
    //        const auto dy = xy[2 * i + 1] - xy[1];
    //        auto d2 = dx * dx + dy * dy;
    //        if (d2 > dmax) {
    //            n1 = i; dmax = d2;
    //        }
    //    }
    //
    //    // 所有点重合
    //    if (dmax < 1E-30)return TP_POINT;
    //
    //    const auto ax = xy[2 * n1 + 0] - xy[0];
    //    const auto ay = xy[2 * n1 + 1] - xy[1];
    //    const auto al = sqrt(dmax);
    //
    //    //--- 2. 计算夹角
    //
    //    for (i = 1; i < np; ++i) {
    //        const auto bx = xy[2 * i + 0] - xy[0];
    //        const auto by = xy[2 * i + 1] - xy[1];
    //        const auto bl = sqrt(bx * bx + by * by);
    //        const auto cos_angle = (ax * bx + ay * by) / (al * bl);
    //        if (fabs(cos_angle) < cos_err)
    //            break;
    //    }
    //
    //    // 夹角在容差范围内，共线
    //    if (i == np) {
    //        //vec3 axis_x{ ax / al,ay / al,0 };
    //        //vec3 origin{ xy[0],xy[1],0 };
    //        //auto axis_y = cross(vec3{ 0,0,1 }, axis_x);
    //        vec2 axis_x{ ax / al,ay / al };
    //        vec2 axis_y{ -axis_x.y, axis_x.x }; // =cross([0,0,1], aixs_x)
    //        vec2 origin{ xy[0],xy[1] };
    //
    //        mt3x3_g2l[0] = axis_x.x; mt3x3_g2l[1] = axis_x.y; mt3x3_g2l[2] = -dot(axis_x, origin);
    //        mt3x3_g2l[3] = axis_y.x; mt3x3_g2l[4] = axis_y.y; mt3x3_g2l[5] = -dot(axis_y, origin);
    //        mt3x3_g2l[6] = 0;        mt3x3_g2l[7] = 0;        mt3x3_g2l[8] = 1;
    //
    //        return ZS_COLINEAR; // 
    //    }
    //    else {
    //        return ZS_COPLANER;
    //    }
    //}

    //! @brief 计算样条插值矩阵中的矩阵TS的逆矩阵，TS矩阵值如下：
    //! 
    //!  (  r00 r01 ... r0N  1  x0  y0  z0 )
    //!  (  r10 r11 ... r1N  1  x1  y1  z1 )
    //!  (  ...                            )
    //!  (  rN0 rN1 ... rNN  1  xN  yN  zN )
    //!  (    1   1 ...   1  0   0   0   0 )
    //!  (   x0  x1 ...  xN  0   0   0   0 )
    //!  (   y0  y1 ...  yN  0   0   0   0 )
    //!  (   z0  z1 ...  zN  0   0   0   0 )
    //! 
    //! @param [in]     np       point number
    //! @param [in]     ndim     dimension, i.e. coordinate number per point
    //! @param [in]     coords   point coordinates, e.g. for 3d: [x0,y0,z0, x1,y1,z1, ...]
    //! @param [out]    mat_ts   Row major squared matrix, rank = ]np + ]ndim + 1
    //! @param [in,out] ibuffer  整数缓冲区，长度 = \np+\ndim+1
    //! @param [out]    singular 矩阵是否奇异，0=否，1=是
    static void compute_xps_ts_inv_(int np, int ndim, const double* coords, double* mat_ts, int* ibuffer, int* singular)
    {
        constexpr double eps = 1E-40;
        const int rank = np + 1 + ndim;// size of matrix TS

        double dx, rij2, * TSi;
        const double* p, * q;
        int i, j, k;

        // 1) left-top block: r^2 log(r^2)
        TSi = mat_ts;
        for (i = 0, p = coords; i < np; ++i) {
            for (j = 0, q = coords; j < np; ++j) {
                // distance
                rij2 = 0;
                for (k = 0; k < ndim; ++k) {
                    dx = p[k] - q[k];
                    rij2 += dx * dx;
                }
                TSi[j] = rij2 * std::log(rij2 + eps);// prevent float overflow

                q += ndim;
            }
            TSi += rank;
            p += ndim;
        }
        // 2) right top block: 1, x, y, z
        for (i = 0, p = coords; i < np; ++i) {
            mat_ts[i * rank + np + 0] = 1;
            for (j = 0; j < ndim; ++j)
                mat_ts[i * rank + np + j + 1] = p[j];
            p += ndim;
        }
        // 3) left bottom block: 
        //   1   1   ...
        //   x0  x1  ...
        //   y0  y1  ...
        //   .   .   ...
        for (j = 0; j < np; ++j) {
            mat_ts[(np + 0) * rank + j] = 1;
            for (k = 0; k < ndim; ++k)
                mat_ts[(np + k + 1) * rank + j] = coords[ndim * j + k];
        }
        // 4) right bottom block: zeros
        for (i = np; i < rank; ++i) {
            for (j = np; j < rank; j++)
                mat_ts[i * rank + j] = 0;
        }

        mat_inverse(rank, mat_ts, ibuffer, singular);
    }

    //! @brief 计算样条插值矩阵
    //! @param [in]     np_src     已知点个数
    //! @param [in]     np_des     未知点（待插值点）个数
    //! @param [in]     ndim       坐标分量个数：1=spline, 2=ips, 3=tps
    //! @param [in]     coords_src 已知点坐标数组
    //! @param [in]     coords_des 未知点坐标数组
    //! @param [in]     ts_inv     TS矩阵的逆，由 \compute_xps_ts_inv 计算得到
    //! @param [out]    mat_G      返回插值矩阵，尺寸为 np_des * np_src
    //! @param [in,out] dbuffer    浮点数缓冲区，长度 = \np_des+\ndim+1
    void compute_xps_interp_matrix_(int np_src, int np_des, int ndim, const double* coords_src, const double* coords_des, const double* ts_inv, double* mat_G, double* dbuffer)
    {
        constexpr double eps = 1E-40;
        const int rank = np_src + 1 + ndim;// size of matrix TS

        double dx, rij2, m, * TFi;
        const double* p, * q, *TS;
        int i, j, k;

        TS  = ts_inv;
        TFi = dbuffer;

        if (!TS || !TFi) { ASSERT(0); }

        // TF = 
        //   r00  r01  ... r0N  1  x0  y0  z0
        //   r10  r11  ... r1N  1  x1  y1  z1
        //   ...
        //   rM0  rM1  ... rMN  1  xM  yM  zM
        //

        //-------------------------------------------------------------------------
        // G = Tf.Ts^-1
        //-------------------------------------------------------------------------
        // NOTE: only output left-top block of G

        // The i-th row of G
        for (i = 0, p = coords_des; i < np_des; ++i) {
            // TF(i) = [rfi1 rfi2 ... rfin 1  xfi yfi ...]
            for (j = 0, q = coords_src; j < np_src; ++j) {
                // distance between source and destination point
                rij2 = 0;
                for (k = 0; k < ndim; ++k) {
                    dx = p[k] - q[k];
                    rij2 += dx * dx;
                }
                // rij
                TFi[j] = rij2 * std::log(rij2 + eps);

                q += ndim;
            }

            // last ndim values
            TFi[np_src + 0] = 1;
            for (k = 0; k < ndim; ++k)
                TFi[np_src + k + 1] = p[k];

            // mat multiply
            for (j = 0; j < np_src; ++j, ++mat_G) {
                m = 0;
                // G(i,j) = TFi(:).TS(:,j)
                for (k = 0; k < rank; ++k) {
                    m += TFi[k] * TS[k * rank + j];
                }
                *mat_G = m;
            }

            p += ndim;
        }
    }

    //------------------------------------------------
    // implements of Boundary
    //------------------------------------------------

    Boundary::Boundary(Boundary&& bd)noexcept
        :nodes_(std::move(bd.nodes_)),
        face_nodes_(std::move(bd.face_nodes_)),
        node_faces_(std::move(bd.node_faces_)),
        node_coords_(std::move(bd.node_coords_)),
        face_centroids_(std::move(bd.face_centroids_)),
        face_types_(std::move(bd.face_types_)),
        face_area_(std::move(bd.face_area_)),
        face_normal_(std::move(bd.face_normal_)),
        topo_(bd.topo_),
        shape_(bd.shape_),
        coord_min_(bd.coord_min_),
        coord_max_(bd.coord_max_),
        mesh_changed_(bd.mesh_changed_),
        //fields_(std::move(bd.fields_)),
        kdtree_(std::move(kdtree_)),
        xps_tm_(bd.xps_tm_),
        xps_coords_(std::move(bd.xps_coords_)),
        xps_ts_inv_(std::move(bd.xps_ts_inv_)),
        xps_computed_(bd.xps_computed_),
        ibuffer_(std::move(bd.ibuffer_)),
        dbuffer_(std::move(bd.dbuffer_)),
        local_xps_points_(std::move(bd.local_xps_points_))
    {}
    Boundary& Boundary::operator = (Boundary&& bd)noexcept
    {
        if (&bd != this) {
            nodes_ = std::move(bd.nodes_);
            face_nodes_ = std::move(bd.face_nodes_);
            node_faces_ = std::move(bd.node_faces_);
            node_coords_ = std::move(bd.node_coords_);
            face_centroids_ = std::move(bd.face_centroids_);
            face_types_ = std::move(bd.face_types_);
            face_area_ = std::move(bd.face_area_);
            face_normal_ = std::move(bd.face_normal_);
            topo_ = bd.topo_;
            shape_ = bd.shape_;
            coord_min_ = bd.coord_min_;
            coord_max_ = bd.coord_max_;
            mesh_changed_ = bd.mesh_changed_;
            //fields_ = std::move(bd.fields_);
            kdtree_ = std::move(kdtree_);
            xps_tm_ = bd.xps_tm_;
            xps_coords_ = std::move(bd.xps_coords_);
            xps_ts_inv_ = std::move(bd.xps_ts_inv_);
            xps_computed_ = bd.xps_computed_;
            ibuffer_ = std::move(bd.ibuffer_);
            dbuffer_ = std::move(bd.dbuffer_);
            local_xps_points_ = std::move(bd.local_xps_points_);
        }        
        return *this;
    }

    void Boundary::clear()
    {
        nodes_.clear();
        face_nodes_.clear();
        node_coords_.clear();
        face_centroids_.clear();
        face_types_.clear();
        face_area_.clear();
        face_normal_.clear();

        topo_  = ZT_POINTS;
        shape_ = ZS_GENERAL;

        xps_tm_.identity();
        xps_coords_.clear();
        xps_ts_inv_.clear();
        xps_computed_ = false;

        kdtree_.clear();
        mesh_changed_ = false;
        is_high_order_ = false;

        fields_.clear();
    }

    void Boundary::reserve(int_l max_node, int_l max_face, int_l max_fnodes)
    {
        max_node = std::max(0, max_node);
        max_face = std::max(0, max_face);
        max_fnodes = std::max(0, max_fnodes);

        nodes_.reserve(max_node);
        node_coords_.reserve(max_node);
        face_centroids_.reserve(max_face);
        face_nodes_.reserve(max_face, max_fnodes);
        //node_faces.reserve(max_node, max_fnodes);
        face_area_.reserve(max_face);
        face_normal_.reserve(max_face);
    }

    int Boundary::add_node(double x, double y, double z, int_g global_id/* = -1*/)
    {
        if (global_id < 0)global_id = (int_g)nodes_.size();

        int id = nodes_.add(global_id);
        
        // add new node
        if (id >= node_coords_.size()) {
            mesh_changed_ = true;
            node_coords_.resize(nodes_.size());
            node_coords_[id].assign(x, y, z);
        }
        // node already exits, modify it's coordinates.
        else {
            set_node_coords(id, x, y, z);
        }

        return id;
    }
    int Boundary::add_node(const vec3& coord, int_g global_id/* = -1*/)
    {
        return add_node(coord.x, coord.y, coord.z, global_id);
    }

    int Boundary::add_face(FaceTopo type, int count, const int_l* fnodes, double cx, double cy, double cz)
    {
        if (type != FT_POLYGON && npf[type] != count) {
            error("Boundary::add_face(), node number not agree.");
            return -1;
        }
        if (type == FT_POLYGON) {
            if      (count == 3)type = FT_TRI3;
            else if (count == 4)type = FT_QUAD4;
        }
        // check
        for (int i = 0; i < count; ++i) {
            if (fnodes[i] < 0 || fnodes[i] >= node_num()) {
                error("Boundary::add_face(), node number out of range.");
                return -1;
            }
            for (int j = i + 1; j < count; ++j) {
                if (fnodes[i] == fnodes[j]) {
                    error("Boundary::add_face(), coincide nodes detected.");
                    return -1;
                }
            }
        }
        if (face_num() == 0)topo_ = ft2zt[type];
        if (ft2zt[type] != topo_) {
            error("Boundary::add_face(), mixed dimension element is not allowed.");
            return -1;
        }

        is_high_order_ &= type == FT_BAR3 || type == FT_TRI6 || type == FT_QUAD8;

        int id = face_nodes_.nrow();
        face_nodes_.push_back(count, fnodes);
        face_types_.push_back(type);
        face_centroids_.emplace_back(cx, cy, cz);
        face_area_.emplace_back(0);
        face_normal_.emplace_back(0, 0, 0);
        
        mesh_changed_ = true;

        return id;
    }
    int Boundary::add_face(FaceTopo type, int count, const int_l* fnodes, const vec3& fcent)
    {
        return add_face(type, count, fnodes, fcent.x, fcent.y, fcent.z);
    }
    int Boundary::add_face(FaceTopo type, int count, const int_l* fnodes)
    {
        return add_face(type, count, fnodes, 0, 0, 0);
    }

    void Boundary::set_node_coords(int_l id, double x, double y, double z)
    {
        auto& c = node_coords_.at(id);
        if (x != c.x || y != c.y || z != c.z)mesh_changed_ = true;
        c.assign(x, y, z);
    }
    void Boundary::set_face_cent(int_l face, double cx, double cy, double cz)
    {
        face_centroids_.at(face).assign(cx, cy, cz);
    }
    void Boundary::set_face_area(int_l face, double sx, double sy, double sz)
    {
        auto ds = std::hypot(sx, sy, sz);
        face_area_.at(face) = ds;
        face_normal_.at(face).assign(sx / ds, sy / ds, sz / ds);
    }

    void Boundary::compute_metics(double biased_angle_deg/* = 5*/)
    {
        if (!mesh_changed_)return;

        // coordinate range
        coord_min_.fill(0);
        coord_max_.fill(0);
        if (node_num() > 0) {
            coord_min_ = coord_max_ = node_coords_.front();
            for (auto& c : node_coords_) {
                coord_min_.x = std::min(coord_min_.x, c.x);
                coord_min_.y = std::min(coord_min_.y, c.y);
                coord_min_.z = std::min(coord_min_.z, c.z);
                coord_max_.x = std::max(coord_max_.x, c.x);
                coord_max_.y = std::max(coord_max_.y, c.y);
                coord_max_.z = std::max(coord_max_.z, c.z);
            }
        }

        // compute face centroid, normal and area.
        if (!face_nodes_.empty()) {
            compute_fc(face_num(), node_coords_.data(), (const FaceTopo*)face_types_.data(), face_nodes_, face_centroids_.data());
            compute_farea(face_num(), node_coords_.data(), (const FaceTopo*)face_types_.data(), face_nodes_, face_normal_.data());
            // normalize
            for (int i = 0; i < face_num(); ++i)face_area_[i] = face_normal_[i].normalize();

            // 
            is_high_order_ = std::all_of(face_types_.begin(), face_types_.end(), [](auto ft) {return ft == FT_BAR3 || ft == FT_TRI6 || ft == FT_QUAD8; });

            // create node-faces
            MeshConnectivity::flip(face_nodes_, node_num(), node_faces_);
        }
        else {
            is_high_order_ = false;
        }

        // shape of boundary
        shape_ = compute_topo_3d(node_num(), node_coords_.data()->data(), xps_tm_.data(), biased_angle_deg);

        // create kdtree
        kdtree_.create(node_coords_.data()->data(), node_num(), true);

        mesh_changed_ = false;
    }

    void Boundary::compute_global_xps_matrix()
    {
        if (xps_computed_)return;

        if (node_num() == 0)return;

        if (mesh_changed_)compute_metics();

        xps_coords_.clear();
        xps_ts_inv_.clear();

        const int_l np = node_num();

        if      (shape_ == ZS_POINT) {
            // do nothing
        }
        else if (shape_ == ZS_COLINEAR) {
            xps_coords_.resize(node_num(), 0);
            for (int_l i = 0; i < node_num(); ++i) {
                auto& p = node_coords_.at(i);
                xps_coords_[i] = xps_tm_(0, 0) * p.x + xps_tm_(0, 1) * p.y + xps_tm_(0, 2) * p.z + xps_tm_(0, 3);
            }
            const int ndim = 1;
            int_l rank = np + ndim + 1;
            
            xps_ts_inv_.resize(rank, rank);
            ibuffer_.resize(rank);
            dbuffer_.resize(rank);

            int singular = 0;
            compute_xps_ts_inv_(node_num(), ndim, xps_coords_.data(), xps_ts_inv_.data(), ibuffer_.data(), &singular);
            if (singular)error("XPS matrix singular, this maybe caused by coincide points!");
        }
        else if (shape_ == ZS_COPLANER) {
            xps_coords_.resize(2 * node_num(), 0);
            for (int_l i = 0; i < node_num(); ++i) {
                auto& p = node_coords_.at(i);
                xps_coords_[2 * i + 0] = xps_tm_(0, 0) * p.x + xps_tm_(0, 1) * p.y + xps_tm_(0, 2) * p.z + xps_tm_(0, 3);
                xps_coords_[2 * i + 1] = xps_tm_(1, 0) * p.x + xps_tm_(1, 1) * p.y + xps_tm_(1, 2) * p.z + xps_tm_(1, 3);
            }

            const int ndim = 2;
            int_l rank = np + ndim + 1;

            xps_ts_inv_.resize(rank, rank);
            ibuffer_.resize(rank);
            dbuffer_.resize(rank);

            int singular = 0;
            compute_xps_ts_inv_(node_num(), 1, xps_coords_.data(), xps_ts_inv_.data(), ibuffer_.data(), &singular);
            if (singular)error("XPS matrix singular, this maybe caused by coincide points!");
        }
        else if (shape_ == ZS_GENERAL) {
            const int ndim = 3;
            int_l rank = np + ndim + 1;
            
            xps_ts_inv_.resize(rank, rank);
            ibuffer_.resize(rank);
            dbuffer_.resize(rank);

            int singular = 0;
            compute_xps_ts_inv_(node_num(), 1, node_coords_.data()->data(), xps_ts_inv_.data(), ibuffer_.data(), &singular);
            if (singular)error("XPS matrix singular, this maybe caused by coincide points!");
        }

        xps_computed_ = true;
    }

    void Boundary::compute_global_xps_interp_coeff(const vec3& p, std::span<double> coeff)
    {
        if (node_num() == 0)return;
        if (!xps_computed_)compute_global_xps_matrix();
        if (coeff.size() < node_num())
            error("invalid length of coefficient array!");

        double q[2]{ 0 };
        const int np_src = static_cast<int>(node_num());
        const int np_des = 1;

        if      (shape_ == ZS_POINT) {
            warn("all boundary points are coincide, the coefficients will be set as 1/n!");
            std::fill(coeff.begin(), coeff.begin() + node_num(), 1.0 / node_num());
        }
        else if (shape_ == ZS_COLINEAR) {
            const int ndim = 1;
            const int_l rank = node_num() + ndim + 1;
            if (dbuffer_.size() < rank)dbuffer_.resize(rank);

            q[0] = xps_tm_(0, 0) * p.x + xps_tm_(0, 1) * p.y + xps_tm_(0, 2) * p.z + xps_tm_(0, 3);
            compute_xps_interp_matrix_(np_src, np_des, ndim, xps_coords_.data(), q, xps_ts_inv_.data(), coeff.data(), dbuffer_.data());
        }
        else if (shape_ == ZS_COPLANER) {
            const int ndim = 2;
            const int_l rank = node_num() + ndim + 1;
            if (dbuffer_.size() < rank)dbuffer_.resize(rank);

            q[0] = xps_tm_(0, 0) * p.x + xps_tm_(0, 1) * p.y + xps_tm_(0, 2) * p.z + xps_tm_(0, 3);
            q[1] = xps_tm_(1, 0) * p.x + xps_tm_(1, 1) * p.y + xps_tm_(1, 2) * p.z + xps_tm_(1, 3);
            compute_xps_interp_matrix_(np_src, np_des, ndim, xps_coords_.data(), q, xps_ts_inv_.data(), coeff.data(), dbuffer_.data());
        }
        else if (shape_ == ZS_GENERAL) {
            const int ndim = 3;
            const int_l rank = node_num() + ndim + 1;
            if (dbuffer_.size() < rank)dbuffer_.resize(rank);

            compute_xps_interp_matrix_(np_src, np_des, ndim, node_coords_.data()->data(), p.data(), xps_ts_inv_.data(), coeff.data(), dbuffer_.data());
        }
    }

    void Boundary::compute_local_xps_interp_coeff(const vec3& p, int max_neigh, std::span<int_l> ids, std::span<double> coeff, int& count)
    {
        if (node_num() == 0)return;
        if (!mesh_changed_)compute_metics();

        if (max_neigh <= 0)error("non-positive max neighbor count!");
        max_neigh = std::min(max_neigh, (int)node_num());
        if (ids.size() < max_neigh || coeff.size() < max_neigh)error("length of ids and coefficients. array is too small!");

        // search nearest points
        if (dbuffer_.size() < max_neigh)dbuffer_.resize(max_neigh);
        count = kdtree_.search(p.data(), max_neigh, ids.data(), coeff.data());
        if (count <= 0)error("failed searching nearest points!");

        // coincide point found
        if (coeff.front() <= 1E-16) {
            count = 1;
            coeff.front() = 1;
            return;
        }

        // extract points
        local_xps_points_.resize(count);
        for (int_l i = 0; i < count; ++i) {
            local_xps_points_[i] = node_coords_[ids[i]];
        }

        // compute shape of nearest point set
        TinyMatrix<double, 4, 4> tm;
        auto type = compute_topo_3d(count, local_xps_points_.data()->data(), tm.data());

        // compute interpolation coefficients

        double q[2]{ 0 };
        int singular = 0;
        if      (type == ZS_POINT   ) {
            ids.front() = 0;
            coeff.front() = 1;
            count = 1;
        }
        else if (type == ZS_COLINEAR) {
            const int ndim   = 1;
            const int np_des = 1;
            const int np_src = count;
            int ni = 0, nd = 0;
            compute_xps_ibuffer_size(np_src, ndim, &ni);
            compute_xps_dbuffer_size(np_src, ndim, &nd);
            
            if (ibuffer_.size() < ni)ibuffer_.resize(ni);
            if (dbuffer_.size() < (ni + ndim * np_des))dbuffer_.resize(ni + ndim * np_des);
            
            auto x_src = dbuffer_.data() + nd;
            ASSERT(x_src + ndim * np_des == dbuffer_.data() + dbuffer_.size());

            // extract local coordinates.
            for (int_l i = 0; i < count; ++i) {
                auto& c = local_xps_points_[i];
                x_src[i] = tm(0, 0) * c.x + tm(0, 1) * c.y + tm(0, 2) * c.z + tm(0, 3);
            }

            // transform query point to local CS.
            q[0] = tm(0, 0) * p.x + tm(0, 1) * p.y + tm(0, 2) * p.z + tm(0, 3);

            // compute interpolation matrix
            compute_xps_interp_matrix(np_src, np_des, ndim, x_src, q, coeff.data(), dbuffer_.data(), ibuffer_.data(), &singular);
        }
        else if (type == ZS_COPLANER) {
            const int ndim   = 2;
            const int np_des = 1;
            const int np_src = count;
            int ni = 0, nd = 0;
            compute_xps_ibuffer_size(np_src, ndim, &ni);
            compute_xps_dbuffer_size(np_src, ndim, &nd);

            if (ibuffer_.size() < ni)ibuffer_.resize(ni);
            if (dbuffer_.size() < (ni + ndim * np_des))dbuffer_.resize(ni + ndim * np_des);

            auto x_src = dbuffer_.data() + nd;
            ASSERT(x_src + ndim * np_des == dbuffer_.data() + dbuffer_.size());

            // extract local coordinates.
            for (int_l i = 0; i < count; ++i) {
                auto& c = local_xps_points_[i];
                x_src[ndim * i + 0] = tm(0, 0) * c.x + tm(0, 1) * c.y + tm(0, 2) * c.z + tm(0, 3);
                x_src[ndim * i + 1] = tm(1, 0) * c.x + tm(1, 1) * c.y + tm(1, 2) * c.z + tm(1, 3);
            }

            // transform query point to local CS.
            q[0] = tm(0, 0) * p.x + tm(0, 1) * p.y + tm(0, 2) * p.z + tm(0, 3);
            q[1] = tm(1, 0) * p.x + tm(1, 1) * p.y + tm(1, 2) * p.z + tm(1, 3);

            // compute interpolation matrix
            compute_xps_interp_matrix(np_src, np_des, ndim, x_src, q, coeff.data(), dbuffer_.data(), ibuffer_.data(), &singular);
        }
        else if (type == ZS_GENERAL ) {
            const int ndim   = 3;
            const int np_des = 1;
            const int np_src = count;

            int ni = 0, nd = 0;
            compute_xps_ibuffer_size(np_src, ndim, &ni);
            compute_xps_dbuffer_size(np_src, ndim, &nd);

            if (ibuffer_.size() < ni)ibuffer_.resize(ni);
            if (dbuffer_.size() < (ni))dbuffer_.resize(ni);

            //x extract neighbor coordinates.

            // compute interpolation matrix
            compute_xps_interp_matrix(np_src, np_des, ndim, local_xps_points_.data()->data(), &p.x, coeff.data(), dbuffer_.data(), ibuffer_.data(), &singular);
        }

        if (singular)error("XPS matrix singular, this maybe caused by coincide points!");
    }

    void Boundary::compute_project_interp_coeff(const vec3& p, int_l(&ids)[npf_max], double(&coeff)[npf_max], int& count)
    {
        if (face_nodes_.empty())error("face data not exists!");

        if (mesh_changed_)compute_metics(); // we need kdtree

        int_l id{ invalid_id };
        double d2{ 0 };
        count = kdtree_.search(p.data(), 1, &id, &d2);
        if (count <= 0)error("failed searching nearest point!");

        // coincide point found
        if (d2 <= 1E-16) {
            ids[0] = id;
            coeff[0] = 1;
            count = 1;
            return;
        }

        using tmat_2x3 = TinyMatrix<double, 2, 3>;
        using tmat_3x3 = TinyMatrix<double, 3, 3>;
        using tmat_4x3 = TinyMatrix<double, 4, 3>;
        using tmat_6x3 = TinyMatrix<double, 6, 3>;
        using tmat_8x3 = TinyMatrix<double, 8, 3>;
        using tmat_1x2 = TinyMatrix<double, 1, 2>;
        using tmat_1x3 = TinyMatrix<double, 1, 3>;
        using tmat_1x4 = TinyMatrix<double, 1, 4>;
        using tmat_1x6 = TinyMatrix<double, 1, 6>;
        using tmat_1x8 = TinyMatrix<double, 1, 8>;

        vec3   Xe[npf_max];
        double Xi[2];
        double Fn[npf_max];
        double d_min = std::numeric_limits<double>::max();
        count = 0;

        // loop each face of node
        for (const auto face : node_faces_[id]) {
            auto nodes = face_nodes_[face];

            // extract points
            for (size_t j = 0; j < nodes.size(); ++j) {
                auto node = nodes[j];
                auto& X = node_coords_[node];
                ASSERT(j < npf_max);
                Xe[j].x = X.x;
                Xe[j].y = X.y;
                Xe[j].z = X.z;
            }

            // projection
            vec3 p_proj;
            switch (face_types_[face]) {
            case FT_BAR2:
                ProjectToLINE2(Xe[0], Xe[1], p, Xi[0]);
                ShapeFunctionLINE2(Xi[0], Fn);
                p_proj = dot(tmat_1x2::view(Fn), tmat_2x3::view(&Xe[0].x));
                break;
            case FT_BAR3:
                ProjectToLINE3(Xe[0], Xe[1], Xe[2], p, Xi[0]);
                ShapeFunctionLINE2(Xi[0], Fn);
                p_proj = dot(tmat_1x3::view(Fn), tmat_3x3::view(&Xe[0].x));
                break;
            case FT_TRI3:
                ProjectToTRI3(Xe, p, Xi[0], Xi[1]);
                ShapeFunctionTRI3(Xi[0], Xi[1], Fn);
                p_proj = dot(tmat_1x3::view(Fn), tmat_3x3::view(&Xe[0].x));
                break;
            case FT_TRI6:
                ProjectToTRI6(Xe, p, Xi[0], Xi[1]);
                ShapeFunctionTRI6(Xi[0], Xi[1], Fn);
                p_proj = dot(tmat_1x6::view(Fn), tmat_6x3::view(&Xe[0].x));
                break;
            case FT_QUAD4:
                ProjectToQUAD4(Xe, p, Xi[0], Xi[1]);
                ShapeFunctionQUAD4(Xi[0], Xi[1], Fn);
                p_proj = dot(tmat_1x4::view(Fn), tmat_4x3::view(&Xe[0].x));
                break;
            case FT_QUAD8:
                ProjectToQUAD8(Xe, p, Xi[0], Xi[1]);
                ShapeFunctionQUAD8(Xi[0], Xi[1], Fn);
                p_proj = dot(tmat_1x8::view(Fn), tmat_8x3::view(&Xe[0].x));
                break;
            default:
                error("invalid face type!");
            }

            // select nearest one
            d2 = distance_sq(p, p_proj);
            if (d2 < d_min) {
                d_min = d2;
                count = static_cast<int>(nodes.size());
                std::copy(nodes.begin(), nodes.end(), ids);
                std::copy(Fn, Fn + count, coeff);
            }
        }
    }
    void Boundary::register_field(const FieldInfo& fd)
    {
        for (auto& f : fields_)if (f.info == &fd)return;
        auto & f = fields_.emplace_back(Field{});
        f.info = &fd;
        f.data.resize(fd.location == NodeCentered ? node_num() : face_num(), fd.ncomp, 0);
    }
    Field& Boundary::get_field(const char* field_name)
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)return f;
        error("field not found!");
        return *(Field*)nullptr;
    }
    const Field& Boundary::get_field(const char* field_name)const
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)return f;
        error("field not found!");
        return *(const Field*)nullptr;
    }

    //void Boundary::register_field(const char* field_name, int ncomp, FieldLocation location, FieldIO iotype, const char* units)
    //{
    //    for(auto& f : fields_){
    //        if (f.name == field_name) {
    //            f.ncomp    = ncomp;
    //            f.location = location;
    //            f.iotype   = iotype;
    //            f.units    = units;
    //
    //            if (location == NodeCentered)
    //                f.data.resize(node_num(), ncomp);
    //            else
    //                f.data.resize(face_num(), ncomp);
    //            return;
    //        }
    //    }
    //
    //    auto& f = fields_.emplace_back(Field{});
    //    f.name     = field_name;
    //    f.ncomp    = ncomp;
    //    f.location = location;
    //    f.iotype   = iotype;
    //    f.units    = units;
    //
    //    if (location == NodeCentered)
    //        f.data.resize(node_num(), ncomp);
    //    else
    //        f.data.resize(face_num(), ncomp);
    //}
    //
    //Field& Boundary::get_field(const char* field_name)
    //{
    //    for (auto& f : fields_)
    //        if (f.name == field_name)return f;
    //
    //    error("field not found!");
    //    return *(Field*)nullptr;
    //}
    //
    //const Field& Boundary::get_field(const char* field_name)const
    //{
    //    for (auto& f : fields_)
    //        if (f.name == field_name)return f;
    //
    //    error("field not found!");
    //    return *(const Field*)nullptr;
    //}
}
