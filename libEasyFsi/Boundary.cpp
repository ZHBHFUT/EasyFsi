#include <cmath>
#include <numbers>  // std::numeric_limits<double>::max();
#include <span>
#include <fstream>  // see Boundary::read_gmsh
#include <set>      // see Boundary::create_edges_

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

    //! @brief ��������㼯�ϵļ���������ʽ�;ֲ�����ϵ
    //! @param np        �����
    //! @param xyz       ���������飬�洢��ʽΪ��[x0,y0,z0, x1,y1,z0, ..., xn,yn,zn]
    //! @param mt4x4_g2l ��ԭ����ϵ��任���ֲ�����ϵ�ı任����size=4*4���������ȴ洢.
    //! @param biased_angle_deg �жϹ�������ƫб�Ƕȣ�0��ʾ�ϸ���.
    //! @return �������˽ṹ���ͣ�0=��(���е��غ�)��1=����(colinear)��2=����(coplaner)��3=һ����ά�ֲ�
    //! @note ��������ʽΪһ����ά�ֲ�ʱ����ı任����Ϊ��λ����.
    static ZoneShape compute_topo_3d(int np, const double* xyz, double* mt4x4_g2l, double biased_angle_deg = 5)
    {
        // �㷨���̣�
        // 1. ����ֲ�����ϵԭ�� origin��ѡȡ��һ����
        // 2. ������ԭ����Զ�ĵ���Ϊ x ���ϵĵ� px����Ҫ�������е�
        // 3. ���� origin��px ���� x �� axis_x��
        // 4. ������ԭ�㡢px����������������ĵ� q����Ҫ�������е�
        // 5. ����origin��px��q����ֲ�����ϵ z �� axis_z
        // 6. ����ֲ�����ϵy�᣺axis_y
        // 7. �γ�����任����

        //
        // �任�����ʽ��
        //          ( ax.x  ax.y  ax.z  -dot(ax,o) )
        //  T_g2l = ( ay.x  ay.y  ay.z  -dot(ay,o) )
        //          ( az.x  az.y  az.z  -dot(az,o) )
        //          (    0     0     0     1       )
        // ���У�oΪ�ֲ�����ϵԭ�����꣬ax��ay��azΪ�ֲ�����ϵ������������ԭ����ϵ�еĵ�λ����ʸ��
        // 
        // ���ǣ�ԭ����ϵ���ͨ�����¹�ʽ�任���ֲ�����ϵ�� Local = T_g2l * Global
        //

        using vec3 = Boundary::vec3;
        using real = double;

        static constexpr auto eps = std::numeric_limits<real>::epsilon();
        static constexpr auto rmax = std::numeric_limits<real>::max();
        const auto cos_err = std::abs(std::cos(std::numbers::pi_v<real> *(90 - biased_angle_deg) / 180));

        //--- 0. ��ʼ������Ϊ��λ����
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                mt4x4_g2l[i * 4 + j] = i != j ? 0 : 1;

        if (np <= 1) { return ZS_POINT; }

        std::span<const vec3> pnts(reinterpret_cast<const vec3*>(xyz), np);

        auto& axis_x = vec3::view(mt4x4_g2l, 0); // x��
        auto& axis_y = vec3::view(mt4x4_g2l, 4); // y��
        auto& axis_z = vec3::view(mt4x4_g2l, 8); // z��

        //--- 1. �ֲ�����ϵԭ�㣺ȡ��һ����

        auto& origin = pnts.front();

        //--- 2. ����x�᣺������ԭ����Զ�ĵ���Ϊ�ֲ�����ϵx���ϵĵ�

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
        // ��λ��x�ᣬ���ģֵΪ0��������е��غ�
        if (axis_x.normalize() <= eps) {
            axis_x.assign(1, 0, 0);
            return ZS_POINT;
        }

        //--- 3. ������ԭ�㡢x����Զ�㹹�ɵ�������������ĵ�

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

        //--- 4. ����ֲ�����ϵ��z�᣺��ֱ�� o-px-q
        axis_z = cross(axis_x, q - origin);
        // ��λ��z�ᣬ���ģֵΪ 0 ��������е㹲��
        if (axis_z.normalize() <= eps) {
            // ѡȡһ����x�᲻��ͬ��ʸ����Ϊz��
            axis_z = cross(axis_x, vec3{ 1,0,0 });
            auto szmax = axis_z.norm_sq();
            auto tz = cross(axis_x, vec3{ 0,1,0 });
            auto sz = tz.norm_sq();
            if (sz > szmax) { szmax = sz; axis_z = tz; }
            tz = cross(axis_x, vec3{ 0,1,0 });
            sz = tz.norm_sq();
            if (sz > szmax) { szmax = sz; axis_z = tz; }

            // ����y��
            axis_y = cross(axis_z, axis_x);

            //--- 5. ���¾ֲ�����ϵ�任����ĵ�4��

            mt4x4_g2l[3] = -dot(axis_x, origin);
            mt4x4_g2l[7] = -dot(axis_y, origin);
            mt4x4_g2l[11] = -dot(axis_z, origin);

            return ZS_COLINEAR;
        }
        // �������������Ƿ���
        else {
            // ����y��
            axis_y = cross(axis_z, axis_x);

            // �������е��� x-y ƽ��ļн�
            for (auto& c : pnts) {
                auto p = c - origin;
                p.normalize();

                // �нǳ�������ݲ������ǹ���
                if (std::abs(dot(axis_z, p)) > cos_err) {
                    axis_x.assign(1, 0, 0);
                    axis_y.assign(0, 1, 0);
                    axis_z.assign(0, 0, 1);
                    return ZS_GENERAL;
                }
            }

            // һЩ��������

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

            //--- 5. ���¾ֲ�����ϵ�任����ĵ�4��

            mt4x4_g2l[3] = -dot(axis_x, origin);
            mt4x4_g2l[7] = -dot(axis_y, origin);
            mt4x4_g2l[11] = -dot(axis_z, origin);

            return ZS_COPLANER;
        }
    }

    //! @brief �������ƽ��㼯�ϵļ���������ʽ�;ֲ�����ϵ
    //! @param np        �����
    //! @param xy        ���������飬�洢��ʽΪ��[x0,y0, x1,y1, ..., xn,yn]
    //! @param mt4x4_g2l ��ԭ����ϵ��任���ֲ�����ϵ�ı任����size=3*3���������ȴ洢.
    //! @param biased_angle_deg �жϹ�������ƫб�Ƕȣ�0��ʾ�ϸ���.
    //! @return �������˽ṹ���ͣ�0=��(���е��غ�)��1=����(colinear)��2=����(coplaner)
    //! @note ��������ʽΪһ���ά�ֲ�ʱ����ı任����Ϊ��λ����.
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
    //    //--- 1. �������1�������Զ�㣬��Ϊ���� a
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
    //    // ���е��غ�
    //    if (dmax < 1E-30)return TP_POINT;
    //
    //    const auto ax = xy[2 * n1 + 0] - xy[0];
    //    const auto ay = xy[2 * n1 + 1] - xy[1];
    //    const auto al = sqrt(dmax);
    //
    //    //--- 2. ����н�
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
    //    // �н����ݲΧ�ڣ�����
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

    //! @brief ����������ֵ�����еľ���TS�������TS����ֵ���£�
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
    //! @param [in,out] ibuffer  ���������������� = \np+\ndim+1
    //! @param [out]    singular �����Ƿ����죬0=��1=��
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

    //! @brief ����������ֵ����
    //! @param [in]     np_src     ��֪�����
    //! @param [in]     np_des     δ֪�㣨����ֵ�㣩����
    //! @param [in]     ndim       �������������1=spline, 2=ips, 3=tps
    //! @param [in]     coords_src ��֪����������
    //! @param [in]     coords_des δ֪����������
    //! @param [in]     ts_inv     TS������棬�� \compute_xps_ts_inv ����õ�
    //! @param [out]    mat_G      ���ز�ֵ���󣬳ߴ�Ϊ np_des * np_src
    //! @param [in,out] dbuffer    ������������������ = \np_des+\ndim+1
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
        face_count_(std::move(bd.face_count_)),
        mesh_changed_(bd.mesh_changed_),
        topo_(bd.topo_),
        shape_(bd.shape_),
        edges_(std::move(bd.edges_)),
        coord_min_(bd.coord_min_),
        coord_max_(bd.coord_max_),
        fields_(std::move(bd.fields_)),
        kdtree_(std::move(kdtree_)),
        xps_tm_(bd.xps_tm_),
        xps_coords_(std::move(bd.xps_coords_)),
        xps_ts_inv_(std::move(bd.xps_ts_inv_)),
        xps_computed_(bd.xps_computed_),
        local_xps_points_(std::move(bd.local_xps_points_)),
        ibuffer_(std::move(bd.ibuffer_)),
        dbuffer_(std::move(bd.dbuffer_))
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
            face_count_ = std::move(bd.face_count_);
            mesh_changed_ = bd.mesh_changed_;
            topo_ = bd.topo_;
            shape_ = bd.shape_;
            edges_ = std::move(bd.edges_);
            coord_min_ = bd.coord_min_;
            coord_max_ = bd.coord_max_;
            fields_ = std::move(bd.fields_);
            kdtree_ = std::move(kdtree_);
            xps_tm_ = bd.xps_tm_;
            xps_coords_ = std::move(bd.xps_coords_);
            xps_ts_inv_ = std::move(bd.xps_ts_inv_);
            xps_computed_ = bd.xps_computed_;
            local_xps_points_ = std::move(bd.local_xps_points_);
            ibuffer_ = std::move(bd.ibuffer_);
            dbuffer_ = std::move(bd.dbuffer_);
        }        
        return *this;
    }

    void Boundary::clear()
    {
        name_.clear();
        nodes_.clear();
        face_nodes_.clear();
        node_coords_.clear();
        face_centroids_.clear();
        face_types_.clear();
        face_area_.clear();
        face_normal_.clear();
        std::fill(face_count_.begin(), face_count_.end(), 0);

        mesh_changed_ = false;

        topo_  = ZT_POINTS;
        shape_ = ZS_GENERAL;

        edges_.clear();
        coord_min_.fill(0);
        coord_max_.fill(0);

        fields_.clear();

        kdtree_.clear();

        xps_tm_.identity();
        xps_coords_.clear();
        xps_ts_inv_.clear();
        xps_computed_ = false;

        local_xps_points_.clear();
        ibuffer_.clear();
        dbuffer_.clear();
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
        face_area_.reserve(max_face);
        face_normal_.reserve(max_face);
    }

    void Boundary::set_name(const char* sname)
    {
        name_ = sname;
        name_.erase(0, name_.find_first_not_of("\r\t\n ")); // trim left
        name_.erase(name_.find_last_not_of("\r\t\n ") + 1); // trim right
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

        int id = face_nodes_.nrow();
        face_nodes_.push_back(count, fnodes);
        face_types_.push_back(type);
        face_centroids_.emplace_back(cx, cy, cz);
        face_area_.emplace_back(0);
        face_normal_.emplace_back(0, 0, 0);
        
        mesh_changed_ = true;

        if (type == FT_POLYGON)edges_.clear();

        ++face_count_[type];

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

    void Boundary::create_edges()
    {
        if (topo_ != ZT_SURFACE || !edges_.empty())return;

        edges_.clear();
        if (face_nodes_.empty() || face_count_[FT_POLYGON] == 0)return;

        // add edge to set
        auto add_edge = [](std::set<Edge>& edges, Edge&& e) {
            if (e.n1 < e.n0) {
                std::swap(e.n0, e.n1);
                std::swap(e.f0, e.f1);
            }

            std::set<Edge>::iterator it = edges.find(e);
            if (it == edges.end()) {
                edges.emplace(e);
            }
            else {
                auto& x = const_cast<Edge&>(*it); //? trick
                if (x.f0 == invalid_id)x.f0 = e.f0;
                else if (x.f1 == invalid_id)x.f1 = e.f1;
            }
        };

        std::set<Edge> edges;
        for (int_l i = 0; i < face_num(); ++i) {
            auto ft = face_types_.at(i);
            auto fnodes = face_nodes_[i];
            switch (ft) {
            case FT_TRI3:
                add_edge(edges, Edge{ fnodes[0], fnodes[1], i, invalid_id });
                add_edge(edges, Edge{ fnodes[1], fnodes[2], i, invalid_id });
                add_edge(edges, Edge{ fnodes[2], fnodes[0], i, invalid_id });
                break;
            case FT_TRI6:
                //         2
                //         +
                //       /  \
                //   5 +     + 4
                //   /        \
                // +-----+-----+
                // 0     3     1
                //
                add_edge(edges, Edge{ fnodes[0], fnodes[3], i, invalid_id });
                add_edge(edges, Edge{ fnodes[3], fnodes[1], i, invalid_id });
                add_edge(edges, Edge{ fnodes[1], fnodes[4], i, invalid_id });
                add_edge(edges, Edge{ fnodes[4], fnodes[2], i, invalid_id });
                add_edge(edges, Edge{ fnodes[2], fnodes[5], i, invalid_id });
                add_edge(edges, Edge{ fnodes[5], fnodes[0], i, invalid_id });
                break;
            case FT_QUAD4:
                add_edge(edges, Edge{ fnodes[0], fnodes[1], i, invalid_id });
                add_edge(edges, Edge{ fnodes[1], fnodes[2], i, invalid_id });
                add_edge(edges, Edge{ fnodes[2], fnodes[3], i, invalid_id });
                add_edge(edges, Edge{ fnodes[3], fnodes[0], i, invalid_id });
                break;
            case FT_QUAD8:
                //        6
                // 3+-----+-----+ 2
                //  |           |
                //  |           |
                // 7+           + 5
                //  |           |
                //  |           |
                //  +-----+-----+
                //  0     4     1
                //
                add_edge(edges, Edge{ fnodes[0], fnodes[4], i, invalid_id });
                add_edge(edges, Edge{ fnodes[4], fnodes[1], i, invalid_id });
                add_edge(edges, Edge{ fnodes[1], fnodes[5], i, invalid_id });
                add_edge(edges, Edge{ fnodes[5], fnodes[2], i, invalid_id });
                add_edge(edges, Edge{ fnodes[2], fnodes[6], i, invalid_id });
                add_edge(edges, Edge{ fnodes[6], fnodes[3], i, invalid_id });
                add_edge(edges, Edge{ fnodes[3], fnodes[7], i, invalid_id });
                add_edge(edges, Edge{ fnodes[7], fnodes[0], i, invalid_id });
                break;
            case FT_POLYGON:
                for (size_t j = 1; j < fnodes.size(); ++j) {
                    add_edge(edges, Edge{ fnodes[j - 1], fnodes[j], i, invalid_id });
                }
                add_edge(edges, Edge{ fnodes.back(), fnodes.front(), i, invalid_id });
                break;
            }
        }

        edges_.reserve(edges.size());
        for (auto& e : edges) {
            ASSERT(e.n0 != e.n1 && e.f0 != e.f1);
            edges_.push_back(e);
        }
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

            // create node-faces
            MeshConnectivity::flip(face_nodes_, node_num(), node_faces_);
        }

        // shape of boundary
        shape_ = compute_topo_3d(node_num(), node_coords_.data()->data(), xps_tm_.data(), biased_angle_deg);

        // create kdtree
        kdtree_.create(node_coords_.data()->data(), node_num(), true);

        // create edges
        if (face_count_[FT_POLYGON] && edges_.empty())create_edges();

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

    void Boundary::compute_local_xps_interp_coeff(const vec3& p, int max_donor, std::span<int_l> ids, std::span<double> coeff, int& count, double min_dist_sq/* = 1E-20*/)
    {
        constexpr const int stride = 3;
        constexpr const int np_des = 1;

        if (node_num() == 0)return;
        if (!mesh_changed_)compute_metics();

        if (max_donor <= 0)error("maximum number of donors is non-positive!");
        max_donor = std::min(max_donor, (int)node_num());
        if (ids.size() < max_donor || coeff.size() < max_donor)error("invalid length of ids and coefficients!");

        // search nearest points
        if (dbuffer_.size() < max_donor)dbuffer_.resize(max_donor);
        count = kdtree_.search(p.data(), max_donor, ids.data(), coeff.data());
        if (count <= 0)error("failed searching nearest points!");

        // coincident point found:
        if (coeff.front() == 0 || coeff.front() <= min_dist_sq) {
            count = 1;
            coeff.front() = 1;
            return;
        }

        const int np_src = count;

        // extract coordinates of donors
        local_xps_points_.resize(np_src);
        for (int_l i = 0; i < np_src; ++i) {
            local_xps_points_[i] = node_coords_[ids[i]];
        }

        // compute shape of donors
        TinyMatrix<double, 4, 4> tm;
        auto type = compute_topo_3d(np_src, local_xps_points_.data()->data(), tm.data());

        int ndim;
        switch (type) {
        case ZS_POINT: // all donors coincided
            ids.front() = 0;
            coeff.front() = 1;
            count = 1;
            return;
        case ZS_COLINEAR:ndim = 1; break;// all donors are colinear
        case ZS_COPLANER:ndim = 2; break;// all donors are coplaner
        default:         ndim = 3;       // general 3-d
        }
        
        // allocate buffer for XPS
        int ni = 0, nd = 0;
        compute_xps_ibuffer_size(np_src, ndim, &ni);
        compute_xps_dbuffer_size(np_src, ndim, &nd);
        if (ibuffer_.size() < ni)ibuffer_.resize(ni);
        if (dbuffer_.size() < nd)dbuffer_.resize(nd);

        // transform donors coordinates to local C.S.
        if (ndim != 3) {
            for (int_l i = 0; i < count; ++i) {
                vec3 c = local_xps_points_[i];
                local_xps_points_[i].x = tm(0, 0) * c.x + tm(0, 1) * c.y + tm(0, 2) * c.z + tm(0, 3);
                local_xps_points_[i].y = tm(1, 0) * c.x + tm(1, 1) * c.y + tm(1, 2) * c.z + tm(1, 3);
                local_xps_points_[i].z = tm(2, 0) * c.x + tm(2, 1) * c.y + tm(2, 2) * c.z + tm(2, 3);
            }
        }

        // transform query point to local C.S.
        vec3 q;
        if (ndim != 3) {
            q.x = tm(0, 0) * p.x + tm(0, 1) * p.y + tm(0, 2) * p.z + tm(0, 3);
            q.y = tm(1, 0) * p.x + tm(1, 1) * p.y + tm(1, 2) * p.z + tm(1, 3);
            q.z = tm(2, 0) * p.x + tm(2, 1) * p.y + tm(2, 2) * p.z + tm(2, 3);
        }
        else
            q = p;

        // compute interpolation matrix
        int singular = 0;
        compute_xps_interp_matrix(np_src, np_des, ndim, stride, local_xps_points_.data()->data(), q.data(), coeff.data(), dbuffer_.data(), ibuffer_.data(), &singular);
        if (!singular)return;

        // singular: using inverse distance method
        error("XPS matrix singular, this is usually caused by coincident donor points!");
        
    }

    void Boundary::compute_project_interp_coeff(const vec3& p, int_l(&ids)[npf_max], double(&coeff)[npf_max], int& count, double min_dist_sq/* = 1E-20*/)
    {
        if (face_nodes_.empty())error("face data not exists!");

        if (mesh_changed_)compute_metics(); // we need kdtree

        int_l id{ invalid_id };
        double d2{ 0 };
        count = kdtree_.search(p.data(), 1, &id, &d2);
        if (count <= 0)error("failed searching nearest point!");

        // coincide point found
        if (d2 ==0 || d2 <= min_dist_sq) {
            ids  [0] = id;
            coeff[0] = 1;
            count    = 1;
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

            if (face_types_[face] == FT_POLYGON) {
                error("Polygon face is unsupported by projection algorithm!");
                return;
            }

            // extract coordinates of donors
            for (size_t j = 0; j < nodes.size(); ++j) {
                auto node = nodes[j];
                auto& X = node_coords_[node];
                ASSERT(j < npf_max);
                Xe[j] = X;
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

            // distance between query point and it's projection
            d2 = distance_sq(p, p_proj);

            // select nearest one
            if (d2 < d_min) {
                d_min = d2;
                count = static_cast<int>(nodes.size());
                std::copy(nodes.begin(), nodes.end(), ids);
                std::copy(Fn, Fn + count, coeff);
            }

            if (d_min == 0 || d_min <= min_dist_sq)break;
        }
    }
    
    void Boundary::register_field_(const FieldInfo& fd)
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

    void Boundary::read_gmsh(const char* file)
    {
        clear();

        std::ifstream fin(file);
        if (!fin) {
            error("%s(), failed open file: %s", __func__, file);
            return;
        }

        enum SectionType {
            ST_Null,
            ST_Comment,
            ST_MeshFormat,
            ST_Nodes,
            ST_Elements,
            ST_Unknown
            //ST_PhysicalNames
        };

        struct ElementData
        {
            int id, type, num_tag, phys_grp, elem_grp, nodes[8];
        };
        //-----------------------------   bar tri quad tet hex prism pyramid
        static const int npe_gmsh[] = { 0, 2, 3,  4,   4,  8,  6,    5 };

        SectionType stype = ST_Null;
        std::string s, s_end;

        while (!fin.eof()) {
            std::getline(fin, s);

            s.erase(0, s.find_first_not_of("\r\t\n "));//�Ƴ���ͷ�Ŀհ�
            s.erase(s.find_last_not_of("\r\t\n ") + 1);//�Ƴ���β�Ŀհ�
            if (s.empty() || s.front() != '$')continue;

            if      (s == "$Comments") {
                stype = ST_Comment;
                s_end = "$EndComments";
            }
            else if (s == "$MeshFormat") {
                stype = ST_MeshFormat;
                s_end = "$EndMeshFormat";
                double version = 0;
                int ftype, size;
                fin >> version >> ftype >> size;
                if (version < 2 || version >= 3) {
                    error("%s(), unsupported GMSH file version: %lf", __func__, version);
                    return;
                }
                if (ftype != 0) {
                    error("%s(), GMSH file is not in ASCII format", __func__);
                    return;
                }
            }
            else if (s == "$Nodes") {
                stype = ST_Nodes;
                s_end = "$EndNodes";

                // NN
                // ID X Y Z
                int nn = 0;
                fin >> nn;
                for (int i = 0; i < nn; ++i) {
                    int id; double x, y, z;
                    fin >> id >> x >> y >> z;
                    add_node(x, y, z, id);
                }
            }
            else if (s == "$Elements") {
                stype = ST_Elements;
                s_end = "$EndElements";

                // Types: bar(1) tri(2) quad(3) tet(4) hex(5) prism(6) pyramid(7)
                //
                // NE
                // Id Type NumTags PhysGrp ElemGrp IndexList

                int nelem = 0;
                fin >> nelem;
                std::vector<ElementData> edata(nelem);
                for (int i = 0; i < nelem; ++i) {
                    int id, type, num_tag, phys_grp, elem_grp, nodes[8];
                    fin >> id >> type >> num_tag >> phys_grp >> elem_grp;
                    if (type < 1 || type >7) {
                        error("%s(), unsupported element type: %d", __func__, type);
                        return;
                    }

                    auto nn = npe_gmsh[type];
                    for (int j = 0; j < nn; ++j)fin >> nodes[j];

                    if (type == 2) {
                        nodes[0] = nodes_.find(nodes[0]).second;
                        nodes[1] = nodes_.find(nodes[1]).second;
                        nodes[2] = nodes_.find(nodes[2]).second;
                        add_face(FT_TRI3, 3, nodes);
                    }
                    else if (type == 3) {
                        nodes[0] = nodes_.find(nodes[0]).second;
                        nodes[1] = nodes_.find(nodes[1]).second;
                        nodes[2] = nodes_.find(nodes[2]).second;
                        nodes[3] = nodes_.find(nodes[3]).second;
                        add_face(FT_QUAD4, 4, nodes);
                    }
                }
            }
            else if (s == s_end) {
                stype = ST_Null;
            }
            else {
                stype = ST_Unknown;
                s_end = "#End" + s.substr(1);
            }
        }
        fin.close();

        mesh_changed_ = true;
        compute_metics();
    }

    std::istream& operator >>(std::istream& is, Boundary& bd)
    {
        bd.clear();

        // nn nf
        int_l nn = 0, nf = 0;
        is >> nn >> nf;
        if (nn <= 0 || nf < 0) {
            error("invalid node or face number in stream!");
            return is;
        }

        bd.reserve(nn, nf, 8 * nf);

        // nodes
        //  ID X Y Z
        for (int_l i = 0; i < nn; ++i) {
            int_g id = 0;
            Boundary::vec3 coords;
            is >> id >> coords.x >> coords.y >> coords.z;
            bd.add_node(coords, id);
        }

        // add faces
        //  TYPE  N1  N2 ...
        std::vector<int_l> fn_l;
        for (int_l i = 0; i < nf; ++i) {
            char type[3] = { '\0' };
            is >> type[0]>>type[1];

            FaceTopo ft = FT_POLYGON;
            if      (_stricmp(type, "L2") == 0)ft = FT_BAR2;
            else if (_stricmp(type, "L3") == 0)ft = FT_BAR3;
            else if (_stricmp(type, "T3") == 0)ft = FT_TRI3;
            else if (_stricmp(type, "T6") == 0)ft = FT_TRI6;
            else if (_stricmp(type, "Q4") == 0)ft = FT_QUAD4;
            else if (_stricmp(type, "Q8") == 0)ft = FT_QUAD8;
            else if (_stricmp(type, "PN") == 0)ft = FT_POLYGON;

            int n2f = 0;
            if (ft != FT_POLYGON)n2f = npf[ft];
            else is >> n2f;

            if (n2f <= 0) {
                error("invalid node number of face in stream!");
                return is;
            }

            fn_l.resize(n2f);
            for (auto& id : fn_l) {
                int_g g = 0;
                is >> g;
                id = bd.nodes().g2l(g);
                if (id == invalid_id) {
                    error("node index out of range in stream!");
                    return is;
                }
            }

            bd.add_face(ft, n2f, fn_l.data());
        }

        // update metrics
        bd.compute_metics();

        return is;
    }
    std::ostream& operator <<(std::ostream& os, const Boundary& bd)
    {
        // nn nf
        os << bd.node_num() << ' ' << bd.face_num() << '\n';

        // nodes
        for (int_l i = 0; i < bd.node_num(); ++i) {
            auto& c = bd.node_coords().at(i);
            os << bd.nodes().l2g(i) << ' '
                << c.x << ' '
                << c.y << ' '
                << c.z << '\n';
        }

        // faces
        for (int_l i = 0; i < bd.face_num(); ++i) {
            auto nodes = bd.face_nodes()[i];
            switch (bd.face_types().at(i)) {
            case FT_BAR2 : os << "L2"; break;
            case FT_BAR3 : os << "L3"; break;
            case FT_TRI3 : os << "T3"; break;
            case FT_TRI6 : os << "T6"; break;
            case FT_QUAD4: os << "Q4"; break;
            case FT_QUAD8: os << "Q8"; break;
            default:
                os << "PN " << nodes.size();
            }

            for (auto id : nodes)os << ' ' << bd.nodes().l2g(id);
            os << '\n';
        }

        return os;
    }
}