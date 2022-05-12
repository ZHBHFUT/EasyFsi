#pragma once
#include "TinyVector.hpp"
#include "TinyMatrix.hpp"

template<typename real = double>
void ShapeFunctionLINE2(real s, real N[2])
{
    N[0] = 0.5 * (1.0 - s);
    N[1] = 0.5 * (1.0 + s);
}
template<typename real = double>
void ShapeFunctionDiff1LINE2(real /*s*/, real dNds[2])
{
    dNds[0] = -1;
    dNds[1] =  1;
}
template<typename real = double>
void ShapeFunctionLINE3(real s, real N[3])
{
    N[0] = (s * s - s) * 0.5;
    N[1] = (s * s + s) * 0.5;
    N[2] = 1.0 - s * s; 
}
template<typename real = double>
void ShapeFunctionDiff1LINE3(real s, real dNds[3])
{
    dNds[0] = s - 0.5;
    dNds[1] = s + 0.5;
    dNds[2] = -2.0 * s;
}
template<typename real = double>
void ShapeFunctionTRI3(real l1, real l2, real N[3])
{
    N[0] = l1;
    N[1] = l2;
    N[2] = 1.0 - l1 - l2;
}
template<typename real = double>
void ShapeFunctionDiff1TRI3(real /*l1*/, real /*l2*/, TinyMatrix<real, 2, 3>& dN)
{
    dN(0, 0) = 1; dN(0, 1) = 0; dN(0, 2) = -1;
    dN(1, 0) = 0; dN(1, 1) = 1; dN(1, 2) = -1;
}
template<typename real = double>
void ShapeFunctionTRI6(real l1, real l2, real N[6])
{
    real l3 = 1 - l1 - l2;
    N[0] = l1 * (2 * l1 - 1);
    N[1] = l2 * (2 * l2 - 1);
    N[2] = l3 * (2 * l3 - 1);
    N[3] = 4 * l1 * l2;
    N[4] = 4 * l2 * l3;
    N[5] = 4 * l3 * l1;
}
template<typename real = double>
void ShapeFunctionDiff1TRI6(real l1, real l2, TinyMatrix<real, 2, 6>& dN)
{
    real l3 = 1 - l1 - l2;

    //         4 * L1 - 1,            0,            1 - 4 * L3,            4 * L2,            -4 * L2,            4 * (L3 - L1)
    dN(0, 0) = 4 * l1 - 1; dN(0, 1) = 0; dN(0, 2) = 1 - 4 * l3; dN(0, 3) = 4 * l2; dN(0, 4) = -4 * l2; dN(0, 5) = 4 * (l3 - l1);
    //         0,            4 * L2 - 1,            1 - 4 * L3,            4 * L1,            4 * (L3 - L2),            -4 * L1
    dN(1, 0) = 0; dN(1, 1) = 4 * l2 - 1; dN(1, 2) = 1 - 4 * l3; dN(1, 3) = 4 * l1; dN(1, 4) = 4 * (l3 - l2); dN(1, 5) = -4 * l1;
}
template<typename real = double>
void ShapeFunctionQUAD4(real s, real t, real N[4])
{
    N[0] = 0.25 * (1.0 - s) * (1.0 - t);
    N[1] = 0.25 * (1.0 + s) * (1.0 - t);
    N[2] = 0.25 * (1.0 + s) * (1.0 + t);
    N[3] = 0.25 * (1.0 - s) * (1.0 + t);
}
template<typename real = double>
void ShapeFunctionDiff1QUAD4(real s, real t, TinyMatrix<real, 2, 4>& dN)
{
    //0.25*[t-1, 1-t, 1+t, -1-t; s-1, -1-s, 1+s, 1-s];
    dN(0, 0) = 0.25 * (t - 1); dN(0, 1) = 0.25 * ( 1 - t); dN(0, 2) = 0.25 * (1 + t); dN(0, 3) = 0.25 * (-1 - t);
    dN(1, 0) = 0.25 * (s - 1); dN(1, 1) = 0.25 * (-1 - s); dN(1, 2) = 0.25 * (1 + s); dN(1, 3) = 0.25 * ( 1 - s);
}
template<typename real = double>
void ShapeFunctionQUAD8(real s, real t, real N[8])
{
    N[0] = 0.25 * (1.0 - s) * (1.0 - t) * (-s - t - 1.0);
    N[1] = 0.25 * (1.0 + s) * (1.0 - t) * ( s - t - 1.0);
    N[2] = 0.25 * (1.0 + s) * (1.0 + t) * ( s + t - 1.0);
    N[3] = 0.25 * (1.0 - s) * (1.0 + t) * (-s + t - 1.0);
    N[4] = 0.5 * (1.0 - s * s) * (1.0 - t);
    N[5] = 0.5 * (1.0 - t * t) * (1.0 + s);
    N[6] = 0.5 * (1.0 - s * s) * (1.0 + t);
    N[7] = 0.5 * (1.0 - t * t) * (1.0 - s);
}
template<typename real = double>
void ShapeFunctionDiff1QUAD8(real s, real t, TinyMatrix<real, 2, 8>& dN)
{
    //0.25*(1-t)*(2*s+t),0.25*(1-t)*(2*s-t),0.25*(1+t)*(2*s+t),0.25*(1+t)*(2*s-t),s*(t-1),0.5*(1-t*t),-s*(t+1),0.5*(t*t-1)
    dN(0, 0) = 0.25 * (1.0 - t) * (2.0 * s + t);
    dN(0, 1) = 0.25 * (1.0 - t) * (2.0 * s - t);
    dN(0, 2) = 0.25 * (1.0 + t) * (2.0 * s + t);
    dN(0, 3) = 0.25 * (1.0 + t) * (2.0 * s - t);
    dN(0, 4) =  s * (t - 1.0);
    dN(0, 5) = 0.5 * (1.0 - t * t);
    dN(0, 6) = -s * (t + 1.0);
    dN(0, 7) = 0.5 * (t * t - 1.0);

    //0.25*(1-s)*(2*t+s),0.25*(1+s)*(2*t-s),0.25*(1+s)*(2*t+s),0.25*(1-s)*(2*t-s),0.5*(s*s-1),-t*(1+s),0.5*(1-s*s),t*(s-1)]
    dN(1, 0) = 0.25 * (1.0 - s) * (2.0 * t + s);
    dN(1, 1) = 0.25 * (1.0 + s) * (2.0 * t - s);
    dN(1, 2) = 0.25 * (1.0 + s) * (2.0 * t + s);
    dN(1, 3) = 0.25 * (1.0 - s) * (2.0 * t - s);
    dN(1, 4) = 0.5 * (s * s - 1.0);
    dN(1, 5) = -t * (1.0 + s);
    dN(1, 6) = 0.5 * (1.0 - s * s);
    dN(1, 7) = t * (s - 1.0);
}

//! project point to line segment.
//! @param [in]  a  start vertex of line
//! @param [in]  b  end vertex of line
//! @param [in]  p  point need to project
//! @param [out] t  parameter, in [-1,1]
//! @note The projection point = (1-s)/2*a + (1+s)/2*b
template<typename real = double>
bool ProjectToLINE2(const TinyVector<real, 3>& a, const TinyVector<real, 3>& b, const TinyVector<real, 3>& p, real& s)
{
    auto ab  = b - a;
    auto dl2 = dot(ab, ab);
    if (dl2 > 0) {
        s = 2.0 * dot(ab, p - a) / dl2 - 1;
        if      (s < -1.0) { s = -1.0; return false; }
        else if (s >  1.0) { s =  1.0; return false; }
        return true;
    }
    else {
        s = -1;
        return false;
    }
}

//! project point to line segment.
//! @param [in]  a  start vertex of line
//! @param [in]  b  end vertex of line
//! @param [in]  c  middle vertex of line
//! @param [in]  p  point need to project
//! @param [out] s  parameter, in [-1,1]
//! @note The projection point = (s*s-s)/2*a + (s*s+s)/2*b + (1-s*s)*c
template<typename real = double>
bool ProjectToLINE3(const TinyVector<real, 3>& a, const TinyVector<real, 3>& b, const TinyVector<real, 3>& c, const TinyVector<real, 3>& p, real& s)
{
    // Q  = Q(s)  --> 线段上的投影点
    // PQ = Q-P   --> 点到投影点的距离矢量
    // d2 = dot(PQ,PQ)
    // d(d2)/ds = 2*dot(d(Q-P)/ds, Q-P)
    //        = 2*dot(dQds, Q-P)
    //        = 2*dot(N'*xyz, N*xyz-P)
    //        = 2*dot(T,PQ)
    //        = 2*f(s) = 0
    // df/ds = dot(N''*xyz, N*xyz-P) + dot(N'*xyz,N'*xyz)
    //       = dot(B, PQ) + dot(T,T)
    //
    // T = N'*xyz
    // B = N''*xyz --> constant
    //
    // N   = [(s^2-s)/2,(s^2+s)/2,1-s^2]
    // N'  = [s-0.5, s+0.5, -2*s]
    // N'' = [1, 1, -2]
    //

    static const size_t NN = 3;//
    static const size_t NC = 3;
    static const size_t NX = 1;

    static const int MAX_ITER = 20;
    static const real ABS_TOL = 1E-20;
    static const real REL_TOL = 1E-10;

    //const real dNds2[3] = { 1,1,-2 };
    //B = dNds2 . xyz;
    TinyVector<real, NN> B {
        1 * a.x + 1 * b.x - 2 * c.x,
        1 * a.y + 1 * b.y - 2 * c.y,
        1 * a.z + 1 * b.z - 2 * c.z
    };
    
    TinyVector<real, NN> ue, Q, T, PQ;
    real f, dfds, ds = 1, dsmax = 0, eps = 1E-8;
    
    s = 0;
    for (int iter = 0; iter < MAX_ITER && std::abs(ds) > eps; ++iter) {
        ShapeFunctionLINE3(s, ue.data());
        Q.x = ue.x * a.x + ue.y * b.x + ue.z * c.x;//ue.xyz
        Q.y = ue.x * a.y + ue.y * b.y + ue.z * c.y;//ue.xyz
        Q.z = ue.x * a.z + ue.y * b.z + ue.z * c.z;//ue.xyz

        ShapeFunctionDiff1LINE3(s, ue.data());
        T.x = ue.x * a.x + ue.y * b.x + ue.z * c.x;//T.xyz
        T.y = ue.x * a.y + ue.y * b.y + ue.z * c.y;//T.xyz
        T.z = ue.x * a.z + ue.y * b.z + ue.z * c.z;//T.xyz

        PQ = Q - p;

        f    = dot(T, PQ);
        dfds = dot(B, PQ) + dot(T, T);
        ds   = -f / dfds;

        dsmax = std::max(std::abs(ds), dsmax);
        eps   = std::max(REL_TOL * dsmax, ABS_TOL);

        s = s + ds;
    }

    //fix
    if      (s < -1.0) { s = -1.0; return false; }
    else if (s >  1.0) { s =  1.0; return false; }
    return true;
}

//! project point to line segment.
//! @param [in]  xyz  3*3 vertex matrix, one row for a vertex
//! @param [in]  p    point need to project
//! @param [out] l1   area coordinate, in [0,1]
//! @param [out] l2   area coordinate, in [0,1]
//! @note The projection point = l1*xyz(0,:) + l2*xyz(1,:) + (1-l1-l2)*xyz(2,:)
template<typename real = double>
bool ProjectToTRI3(const TinyVector<real, 3> Xe[3], const TinyVector<real, 3>& p, real& l1, real& l2)
{
    constexpr real eps = std::numeric_limits<real>::epsilon();

    //
    //           V2
    //           +
    //          / \
    //        /    \
    //      /       \
    //     +---------+
    //    V0         V1
    //

    static const size_t NN = 3;//
    static const size_t NC = 3;
    static const size_t NX = 2;

    static const TinyVector<real, NN> v0{ 0,0,0 };

    auto v1 = Xe[1] - Xe[0];
    auto v2 = Xe[2] - Xe[0];
    auto vp = p     - Xe[0];

    real y22z22   = v1.y*v1.y + v1.z*v1.z;
    real y32z32   = v2.y*v2.y + v2.z*v2.z;
    real y3z2y2z3 = v2.y*v1.z - v1.y*v2.z;
    real y2y3z2z3 = v1.y*v2.y + v1.z*v2.z;

    real bb = v2.x*v2.x * y22z22 + y3z2y2z3*y3z2y2z3 - 2 * v1.x*v2.x*y2y3z2z3 + v1.x*v1.x*y32z32;
    
    real l3 = 0;
    if (bb < -eps || bb > eps) {
        // !!!NOTE THIS!!!
        l2 = (v2.x * v2.x * (vp.y * v1.y + vp.z * v1.z) + (v2.y * vp.z - vp.y * v2.z) * y3z2y2z3 - v2.x * (v1.x * (vp.y * v2.y + vp.z * v2.z) + vp.x * y2y3z2z3) + vp.x * v1.x * y32z32) / bb;
        l3 = (v1.x * v1.x * (vp.y * v2.y + vp.z * v2.z) - (v1.y * vp.z - vp.y * v1.z) * y3z2y2z3 - v1.x * (v2.x * (vp.y * v1.y + vp.z * v1.z) + vp.x * y2y3z2z3) + vp.x * v2.x * y22z22) / bb;
    }
    else {
        // colinear
        ProjectToLINE2(v0, v1, vp, l2);//l2:[-1,1]
        l2 = 0.5 * (l2 + 1.0);//to [0,1]
        if (l2 == -1) {
            ProjectToLINE2(v1, v2, vp, l3);//l3:[-1,1]
            l3 = 0.5 * (l3 + 1.0);//to [0,1]
            l2 = 1 - l3;
            l1 = 0;
        }
        else
            l1 = 1 - l2;
        return false;
    }
    
    l1 = 1 - l2 - l3;
    
    //TRI
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
            //project to line(v0,v1)
            ProjectToLINE2(v0, v1, vp, l2);//l2:[-1,1]
            l2 = 0.5 * (l2 + 1.0);//to [0,1]
            l1 = 1 - l2;
        }
    }
    else if (l1 <= 0) {
        if (l2 <= 0) {
            l1 = l2 = 0;
        }
        else {
            //project to line(v1,v2)
            ProjectToLINE2(v1, v2, vp, l3);//l3:[-1,1]
            l3 = 0.5 * (l3 + 1.0);//to [0,1]
            l2 = 1 - l3;
            l1 = 0;
        }
    }
    else {
        //project to line(v0,v2)
        ProjectToLINE2(v0, v2, vp, l3);//l3:[-1,1]
        l3 = 0.5 * (l3 + 1.0);//to [0,1]
        l1 = 1 - l3;
        l2 = 0;
    }
    return false;
}

//! project point to line segment.
//! @param [in]  Xe   6*3 vertex matrix, one row for a vertex
//! @param [in]  p    point need to project
//! @param [out] l1   area coordinate, in [0,1]
//! @param [out] l2   area coordinate, in [0,1]
//! @note The projection point = l1*xyz(0,:) + l2*xyz(1,:) + (1-l1-l2)*xyz(2,:)
template<typename real = double>
bool ProjectToTRI6(const TinyVector<real, 3> Xe[6], const TinyVector<real, 3>& p, real& l1, real& l2)
{
    //
    //           V2
    //           +
    //          / \
    //     V5 +    + V4
    //      /       \
    //     +----+----+
    //    V0    V3   V1
    //

    static const size_t NN = 6;//
    static const size_t NC = 3;
    static const size_t NX = 2;

    static const int MAX_ITER = 20;
    static const real ABS_TOL = 1E-20;
    static const real REL_TOL = 1E-10;

    static const TinyVector<real, NN> A11{ 4,0,4,0, 0,-8 };//d(dNdL1)/dL1
    static const TinyVector<real, NN> A12{ 0,0,4,4,-4,-4 };//d(dNdL1)/dL2 
    static const TinyVector<real, NN> A21{ 0,0,4,4,-4,-4 };//d(dNdL2)/dL1 = d(dNdL1)/dL2 
    static const TinyVector<real, NN> A22{ 0,4,4,0,-8, 0 };//d(dNdL2)/dL2

    const auto& xyz = TinyMatrix<real, 6, 3>::view(Xe->data());

    auto B11 = dot(A11, xyz);
    auto B12 = dot(A12, xyz);
    auto B21 = dot(A21, xyz);
    auto B22 = dot(A22, xyz);
    
    TinyVector<real, NN> N;
    TinyVector<real, NC> Q, PQ;
    TinyMatrix<real, NX, NN> dN;
    TinyMatrix<real, NX, NC> T;
    TinyVector<real, NX> f, dl;
    TinyMatrix<real, NX> J;
    real l3, dl1max = 0, dl2max = 0, dlsqr, relsqr;

    l1 = 1.0 / 3.0, l2 = 1.0 / 3.0;
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        //l3 = 1.0 - l1 - l2;
        ShapeFunctionTRI6(l1, l2, N.data());
        Q = dot(N, xyz);

        ShapeFunctionDiff1TRI6(l1, l2, dN);
        T = dot(dN, xyz);

        PQ = Q - p;

        f = dot(T, PQ);

        //
        J(0, 0) = dot(B11, PQ); J(0, 1) = dot(B12, PQ);
        J(1, 0) = dot(B21, PQ); J(1, 1) = dot(B22, PQ);

        //J + T.T'
        auto Tt = transpose(T);
        J += dot(T, Tt);

        //dl = -J\f
        inverse(J);
        dl = -dot(J, f);

        l1 += dl[0];
        l2 += dl[1];

        //converge check

        dl1max = std::max(std::abs(dl[0]), dl1max);
        dl2max = std::max(std::abs(dl[1]), dl2max);

        dlsqr  = dot(dl, dl);
        real r0 = dl[0] / dl1max;
        real r1 = dl[1] / dl2max;
        relsqr = r0 * r0 + r1 * r1;

        if (dlsqr  <= ABS_TOL * ABS_TOL ||
            relsqr <= REL_TOL * REL_TOL)break;
    }

    l3 = 1.0 - l1 - l2;

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
            ProjectToLINE3(Xe[0], Xe[1], Xe[3], p, l2);//l2:[-1,1]
            l2 = 0.5 * (l2 + 1.0);//to [0,1]
            l1 = 1 - l2;
        }
    }
    else if (l1 <= 0) {
        if (l2 <= 0) {
            l1 = l2 = 0;
        }
        else {
            //project to line(v1,v2,v4)
            ProjectToLINE3(Xe[1], Xe[2], Xe[4], p, l3);//l3:[-1,1]
            l3 = 0.5 * (l3 + 1.0);//to [0,1]
            l2 = 1 - l3;
            l1 = 0;
        }
    }
    else {
        //project to line(v0,v2,v5)
        ProjectToLINE3(Xe[0], Xe[2], Xe[5], p, l3);//l3:[-1,1]
        l3 = 0.5 * (l3 + 1.0);//to [0,1]
        l1 = 1 - l3;
        l2 = 0;
    }
    return false;
}

//! project point to line segment.
//! @param [in]  Xe  4*3 vertex matrix, one row for a vertex
//! @param [in]  p   point need to project
//! @param [out] s   natural coordinate, in [-1,1]
//! @param [out] t   natural coordinate, in [-1,1]
template<typename real = double>
bool ProjectToQUAD4(const TinyVector<real, 3> Xe[4], const TinyVector<real, 3>& p, real& s, real& t)
{
    //
    //   V3+---------+ V2
    //     |         |
    //     |         |
    //     |         |
    //     |         |
    //     +---------+
    //    V0         V1
    //

    static const size_t NN = 4;//
    static const size_t NC = 3;
    static const size_t NX = 2;

    static const int MAX_ITER = 20;
    static const real ABS_TOL = 1E-20;
    static const real REL_TOL = 1E-10;

    static const TinyVector<real, NN> A11{ 0,0,0,0 };//d(dNds)/ds
    static const TinyVector<real, NN> A12{ 0.25,-0.25,0.25,-0.25 };//d(dNds)/dt 
    static const TinyVector<real, NN> A21{ 0.25,-0.25,0.25,-0.25 };//d(dNdt)/ds = d(dNds)/dt
    static const TinyVector<real, NN> A22{ 0,0,0,0 };//d(dNdt)/dt

    const auto& xyz = TinyMatrix<real, NN, NC>::view(Xe->data());

    auto B11 = dot(A11, xyz);
    auto B12 = dot(A12, xyz);
    auto B21 = dot(A21, xyz);
    auto B22 = dot(A22, xyz);

    TinyVector<real, NN> N;
    TinyVector<real, NC> Q, PQ;
    TinyMatrix<real, NX, NN> dN;
    TinyMatrix<real, NX, NC> T;
    TinyVector<real, NX> f, dX;
    TinyMatrix<real, NX> J;
    
    real dl1max = 0, dl2max = 0, dlsqr, relsqr;

    s = 0; t = 0;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        ShapeFunctionQUAD4(s, t, N.data());
        Q = dot(N, xyz);

        ShapeFunctionDiff1QUAD4(s, t, dN);
        T = dot(dN, xyz);

        PQ = Q - p;

        f = dot(T, PQ);

        //
        J(0, 0) = dot(B11, PQ); J(0, 1) = dot(B12, PQ);
        J(1, 0) = dot(B21, PQ); J(1, 1) = dot(B22, PQ);

        //J + T.T'
        auto Tt = transpose(T);
        J += dot(T, Tt);

        //dl = -J\f
        inverse(J);
        dX = -dot(J, f);

        s += dX[0];
        t += dX[1];

        //converge check

        dl1max = std::max(std::abs(dX[0]), dl1max);
        dl2max = std::max(std::abs(dX[1]), dl2max);

        dlsqr = dot(dX, dX);
        real r0 = dX[0] / dl1max;
        real r1 = dX[1] / dl2max;
        relsqr = r0 * r0 + r1 * r1;

        if (dlsqr <= ABS_TOL * ABS_TOL ||
            relsqr <= REL_TOL * REL_TOL)break;
    }
    

    if (s < -1) {
        ProjectToLINE2(Xe[0], Xe[3], p, t);
        s = -1;
        return false;
    }
    else if (s > 1) {
        ProjectToLINE2(Xe[1], Xe[2], p, t);
        s = 1;
        return false;
    }
    else if (t < -1) {
        ProjectToLINE2(Xe[0], Xe[1], p, s);
        t = -1;
        return false;
    }
    else if (t > 1) {
        ProjectToLINE2(Xe[3], Xe[2], p, s);
        t = 1;
        return false;
    }
    return true;
}

//! project point to line segment.
//! @param [in]  Xe  8*3 vertex matrix, one row for a vertex
//! @param [in]  p   point need to project
//! @param [out] s   natural coordinate, in [-1,1]
//! @param [out] t   natural coordinate, in [-1,1]
template<typename real = double>
bool ProjectToQUAD8(const TinyVector<real, 3> Xe[8], const TinyVector<real, 3>& p, real& s, real& t)
{
    //          V6
    //   V3+----+----+ V2
    //     |         |
    //     |         |
    //  V7 +         + V5
    //     |         |
    //     |         |
    //     +----+----+
    //    V0    V4   V1
    //

    static const size_t NN = 8;//
    static const size_t NC = 3;
    static const size_t NX = 2;

    static const int MAX_ITER = 20;
    static const real ABS_TOL = 1E-20;
    static const real REL_TOL = 1E-10;

    const auto& xyz = TinyMatrix<real, NN, NC>::view(Xe->data());

    //d(dNds)/ds
    auto getA11 = [](real /*s*/, real t, TinyVector<real, NN>& A11) {
        //0.5*(1-t),0.5*(1-t),0.5*(1+t),0.5*(1+t),-1+t,0,-1-t,0
        A11[0] = 0.5 * (1 - t);
        A11[1] = 0.5 * (1 - t);
        A11[2] = 0.5 * (1 + t);
        A11[3] = 0.5 * (1 + t);
        A11[4] = -1.0 + t;
        A11[5] = 0.0;
        A11[6] = -1.0 - t;
        A11[7] = 0.0;
    };
    //d(dNds)/dt
    auto getA12 = [](real s, real t, TinyVector<real, NN>& A12) {
        //0.25*(1-2*s-2*t),0.25*(-1-2*s+2*t),0.25*(1+2*s+2*t),0.25*(-1+2*s-2*t),s,-t,-s,t
        A12[0] = 0.25 * ( 1.0 - 2.0 * s - 2.0 * t);
        A12[1] = 0.25 * (-1.0 - 2.0 * s + 2.0 * t);
        A12[2] = 0.25 * ( 1.0 + 2.0 * s + 2.0 * t);
        A12[3] = 0.25 * (-1.0 + 2.0 * s - 2.0 * t);
        A12[4] =  s;
        A12[5] = -t;
        A12[6] = -s;
        A12[7] =  t;
    };
    //d(dNdt)/dt
    auto getA22 = [](real s, real /*t*/, TinyVector<real, NN>& A22) {
        //0.5*(1-s),0.5*(1+s),0.5*(1+s),0.5*(1-s),0,-1-s,0,-1+s
        A22[0] = 0.5 * (1 - s);
        A22[1] = 0.5 * (1 + s);
        A22[2] = 0.5 * (1 + s);
        A22[3] = 0.5 * (1 - s);
        A22[4] = 0.0;
        A22[5] = -1.0 - s;
        A22[6] = 0.0;
        A22[7] = -1.0 + s;
    };

    TinyVector<real, NN> N;
    TinyVector<real, NC> Q, PQ;
    TinyMatrix<real, NX, NN> dN;
    TinyMatrix<real, NX, NC> T;
    TinyVector<real, NX> f, dX;
    TinyMatrix<real, NX> J;
    TinyVector<real, NC> B11, B12, B22;
    TinyVector<real, NN> A11;//d(dNds)/ds
    TinyVector<real, NN> A12;//d(dNds)/dt 
    //TinyVector<real, NN> A21;//d(dNdt)/ds = d(dNds)/dt
    TinyVector<real, NN> A22;//d(dNdt)/dt
    real dsmax = 0, dtmax = 0, dlsqr, relsqr;

    s = 0; t = 0;
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        ShapeFunctionQUAD8(s, t, N.data());
        Q = dot(N, xyz);

        ShapeFunctionDiff1QUAD8(s, t, dN);
        T = dot(dN, xyz);

        getA11(s, t, A11);
        getA12(s, t, A12);
        getA22(s, t, A22);

        B11 = dot(A11, xyz);
        B12 = dot(A12, xyz);
        //B21 = B12;// dot(A21, xyz);
        B22 = dot(A22, xyz);

        PQ = Q - p;

        f = dot(T, PQ);

        //
        J(0, 0) = dot(B11, PQ); J(0, 1) = dot(B12, PQ);
        J(1, 0) = dot(B12, PQ); J(1, 1) = dot(B22, PQ);

        //J + T.T'
        auto Tt = transpose(T);
        J += dot(T, Tt);

        //dl = -J\f
        inverse(J);
        dX = -dot(J, f);

        s += dX[0];
        t += dX[1];

        //converge check

        dsmax = std::max(std::abs(dX[0]), dsmax);
        dtmax = std::max(std::abs(dX[1]), dtmax);

        dlsqr = dot(dX, dX);
        real r0 = dX[0] / dsmax;
        real r1 = dX[1] / dtmax;
        relsqr = r0 * r0 + r1 * r1;

        if (dlsqr <= ABS_TOL * ABS_TOL ||
            relsqr <= REL_TOL * REL_TOL)break;
    }

    if (s < -1) {
        ProjectToLINE3(Xe[0], Xe[3], Xe[7], p, t);
        s = -1;
        return false;
    }
    else if (s > 1) {
        ProjectToLINE3(Xe[1], Xe[2], Xe[5], p, t);
        s = 1;
        return false;
    }
    else if (t < -1) {
        ProjectToLINE3(Xe[0], Xe[1], Xe[4], p, s);
        t = -1;
        return false;
    }
    else if (t > 1) {
        ProjectToLINE3(Xe[3], Xe[2], Xe[6], p, s);
        t = 1;
        return false;
    }
    return true;
}

/*
//-------------------------------
// Test
//-------------------------------

//QUAD8
{
    real s = 0, t = 0;
    TinyVector<real, 3> Xe[8] = {
        {-1,-1,0},
        { 1,-1,0},
        { 1, 1,0},
        {-1, 1,0},
        { 0,-1,0},
        { 1, 0,0},
        { 0, 1,0},
        {-1, 0,0}
    };
    TinyVector<real, 3> p{ 0.5,0.25,1.0 };
    ProjectToQUAD8(Xe, p, s, t);//s=0.5, t=0.25, return=true
    printf("s = %g, t = %g\n", s, t);
    p = { -2.0,0.5,1.0 };
    ProjectToQUAD8(Xe, p, s, t);//s=-1, t=0.5, return=false
    printf("s = %g, t = %g\n", s, t);
    p = { -2.0,-2.0,1.0 };
    ProjectToQUAD8(Xe, p, s, t);//s=-1, t=-1, return=false
    printf("s = %g, t = %g\n", s, t);
}

*/