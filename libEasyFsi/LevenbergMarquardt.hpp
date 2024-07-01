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
//! \file LevenbergMarquardt.hpp
//! \brief LevenbergMarquardt algorithm for solving non-linear least-squares problems.
//! \author ZHANG Bing, zhangbing@hfut.edu.cn
//! \date 2024-04-18
//! \copyright (c)2024 All rights reserved.
//!-----------------------------------------------------------------------------

#include <limits>

#include "VecMatAlg.hpp"

namespace EasyLib {

    //--------------------------------------------------------
    // Solve non-linear least-squares problems:
    //   min(||x - f(p)||)
    // with Levenberg-Marquardt Algorithm.
    // 
    // where,
    //    x    = {x1,x2,...,xM} is a constant vector.
    //    f(p} = {f1,f2,...,fN} is a multivariate function.
    //    p    = {p1,p2,...,pM} is vector of independent variables.
    //
    // reference: 
    // [1] http://users.ics.forth.gr/~lourakis/levmar/levmar.pdf
    //
    //--------------------------------------------------------

    //! @brief Levenberg-Marquardt Algorithm for solving: min(||x - f(p)||).
    //! @tparam M The dimension of p, i.e. the unknown number.
    //! @tparam N The dimension of f and x. N >= M.
    template<std::size_t M, std::size_t N>
    struct MinimizeLM
    {
        static_assert(M > 0 && N >= M);

        // minimizing: ||x-f(p)||
        template<
            typename ObjectiveFunction = void(const VecMatAlg::tvec<double, M>& p, VecMatAlg::tvec<double, N>& f),
            typename JacobianFunction = void(const VecMatAlg::tvec<double, M>& p, VecMatAlg::tmat<double, N>& J),
            typename VecF = VecMatAlg::tvec<double, N>,
            typename VecP = VecMatAlg::tvec<double, M>
        >
        static void solve(ObjectiveFunction getF, JacobianFunction getJ, const VecF& x, VecP& p0, int max_iter = 100)
        {
            using namespace VecMatAlg;
            using type = std::remove_cv_t<std::remove_reference_t<decltype(p0[0])>>;
            using matJ = tmat<type, N, M>; //? fix me: Jacobian should be N*3
            using vecP = tvec<type, M>;
            using vecF = tvec<type, N>;
            using matH = tmat<type, M>;

            constexpr auto epsilon   = std::numeric_limits<type>::epsilon();
            constexpr type tau       = 1E-3;
            constexpr type one_three = type{ 1 } / type{ 3 };

            auto& vec_p0 = vec_view<type, M>(&p0[0]);
            auto& vec_x  = vec_view<type, N>(& x[0]);

            //
            type nu = 2;
            vecP p  = vec_p0;

            vecF f; // N*1
            getF(p, f);

            // J = {{df1/dp1,df1/dp2,...,df1/dpM},
            //      {df2/dp1,df2/dp2,...,df2/dpM},
            //      ...
            //      {dfN/dp1,dfN/dp2,...,dfN/dpM}}
            matJ J; // N*M
            getJ(p, J);

            matH A = transpose_dot(J, J); // M*M, = J^T . J

            auto eps_p       = vec_x - f; // N*1
            auto eps_p_norm2 = norm_sq(eps_p); // = ||eps_p||^2
            auto g           = transpose_dot(J, eps_p); // = J^T . eps_p

            auto stop = norm_infinity(g) < epsilon;

            auto max_Aii = A[0][0];
            for (std::size_t i = 1; i < M; ++i)max_Aii = std::max(max_Aii, A[i][i]);
            auto mu = tau * max_Aii;

            int it = 0, n = 0;
            for (; !stop && it < max_iter; ++it) {
                for (; !stop;) {
                    ++n;

                    // solve: (A + mu I) dp = g
                    for (std::size_t i = 0; i < M; ++i)A[i][i] += mu;
                    vecP dp{ 0 };
                    if (!cholesky_solve(A, g, dp))
                        lu_solve(A, g, dp);

                    // stop = ||dp|| <= eps2 * ||p||
                    stop = norm_sq(dp) <= epsilon * epsilon * norm_sq(p);

                    if (!stop) {
                        // p_new = p + dp
                        auto p_new = p + dp;

                        // rho = (||eps_p||^2 - ||x-f||^2)/(dp' . (mu * dp + g))
                        getF(p_new, f);
                        auto rho = (eps_p_norm2 - distance_sq(vec_x, f)) / (dot(dp, mu * dp + g));
                        if (rho > 0) {
                            p = p_new;
                            getJ(p, J);

                            A = transpose_dot(J, J); // M*M, = J^T . J

                            eps_p       = vec_x - f;
                            eps_p_norm2 = norm_sq(eps_p); // = ||eps_p||^2

                            g = transpose_dot(J, eps_p); // = J^T . eps_p

                            stop = norm_infinity(g) < epsilon || eps_p_norm2 < epsilon * epsilon;
                            auto t = 2 * rho - 1;
                            mu *= std::max(one_three, 1 - t * t * t);
                            nu = 2;
                            break;
                        }
                        else {
                            mu *= nu;
                            nu *= 2;
                        }
                    }
                }
            }

            vec_p0 = p;
        }
    };

    template<>
    struct MinimizeLM<1, 1>
    {
        //! @brief Find solution of min(||x - f(p)||).
        //! @tparam ObjectiveFunction The type of multivariate function f: void(auto p, auto& f);
        //! @tparam JacobianFunction  The function used to evaluate Jacobian matrix, i.e. df/dp: void(auto x, auto& J);
        //! @tparam VecF  The vector type of multivariate function f.
        //! @tparam VecX  The vector type of independent variables p.
        //! @param [in]    getF  The multivariate function f.
        //! @param [in]    getJ  The function used to evaluate Jacobian matrix, i.e. df/dp
        //! @param [inout] x0    The initial and solution values of independent variables x.
        //! @param [in]    max_iter Maximum iteration number, default is 100.
        template<
            typename ObjectiveFunction = void(double p, double& f),
            typename JacobianFunction  = void(double p, double& J),
            typename VecF = double,
            typename VecP = double
        >
        static void solve(ObjectiveFunction getF, JacobianFunction getJ, const VecF x, VecP& p0, int max_iter = 100)
        {
            using namespace VecMatAlg;
            using type = std::remove_cv_t<std::remove_reference_t<decltype(p0)>>;
            
            constexpr auto epsilon   = std::numeric_limits<type>::epsilon();
            constexpr auto tau       = type{ 1E-3 };
            constexpr auto one_three = type{ 1 } / type{ 3 };

            //
            type nu = 2;
            type p  = p0;

            type J = 0;
            type f = 0;
            getF(p, f);
            getJ(p, J);

            auto A     = J * J;
            auto eps_p = x - f;
            auto eps_p_norm2 = eps_p * eps_p;
            auto g     = J * eps_p;
            auto stop  = std::fabs(g) < epsilon;
            auto mu    = tau * A;

            int it = 0, n = 0;
            for (; !stop && it < max_iter; ++it) {
                for (; !stop;) {
                    ++n;

                    // solve: (A + mu I) dp = g
                    auto dp = g / (A + mu);

                    // stop = ||dp| <= eps2 * ||p||
                    stop = std::fabs(dp) < epsilon * std::fabs(p);

                    if (!stop) {
                        // p_new = p + dp
                        auto p_new = p + dp;

                        // rho = (||eps_p||^2 - ||x-f||^2)/(dp' . (mu * dp + g))
                        getF(p_new, f);
                        auto rho = (eps_p_norm2 - (x - f) * (x - f)) / (dp * (mu * dp + g));
                        if (rho > 0) {
                            p = p_new;
                            getJ(p, J);
                            A = J * J; // = J' . J
                            eps_p = x - f;
                            eps_p_norm2 = eps_p * eps_p; // = ||eps_p||^2
                            g = J * eps_p; // = J' . eps_p

                            stop = std::fabs(g) < epsilon || eps_p_norm2 < epsilon * epsilon;

                            auto t = 2 * rho - 1;
                            mu *= std::max(one_three, 1 - t * t * t);
                            nu = 2;
                            break;
                        }
                        else {
                            mu *= nu;
                            nu *= 2;
                        }
                    }
                }
            }

            p0 = p;
        }
    };

    template<std::size_t N>
    struct MinimizeLM<1, N>
    {
        static_assert(N > 1);

        // minimizing: ||x-f(p)||
        template<
            typename ObjectiveFunction = void(double p, VecMatAlg::tvec<double, N>& f),
            typename JacobianFunction  = void(double p, VecMatAlg::tmat<double, N>& J),
            typename VecF              = VecMatAlg::tvec<double, N>,
            typename VecP              = double
        >
        static void solve(ObjectiveFunction getF, JacobianFunction getJ, const VecF& x, VecP& p0, int max_iter = 100)
        {
            using namespace VecMatAlg;
            using type = std::remove_cv_t<std::remove_reference_t<decltype(p0)>>;
            using matJ = tvec<type, N>;
            using vecF = tvec<type, N>;

            constexpr auto epsilon   = std::numeric_limits<type>::epsilon();
            constexpr auto tau       = type{ 1E-3 };
            constexpr auto one_three = type{ 1 } / type{ 3 };

            auto& vec_x = vec_view<type, N>(&x[0]);

            //
            type nu = 2;
            type p  = p0;

            matJ J;
            vecF f;
            getF(p, f);
            getJ(p, J);

            auto A     = dot(J, J);
            auto eps_p = vec_x - f;
            auto eps_p_norm2 = norm_sq(eps_p);
            auto g    = dot(J, eps_p);
            auto stop = std::fabs(g) < epsilon;
            auto mu   = tau * A;
        
            int it = 0, n = 0;
            for (; !stop && it < max_iter; ++it) {
                for (; !stop;) {
                    ++n;

                    // solve: (A + mu I) dp = g
                    auto dp = g / (A + mu);

                    // stop = ||dp| <= eps2 * ||p||
                    stop = std::fabs(dp) < epsilon * std::fabs(p);

                    if (!stop) {
                        // p_new = p + dp
                        auto p_new = p + dp;

                        // rho = (||eps_p||^2 - ||x-f(p)||^2)/(dp' . (mu * dp + g))
                        getF(p_new, f);
                        auto rho = (eps_p_norm2 - distance_sq(vec_x, f)) / (dp * (mu * dp + g));
                        if (rho > 0) {
                            p = p_new;
                            getJ(p, J);

                            A = dot(J, J); // = J' . J

                            eps_p       = vec_x - f;
                            eps_p_norm2 = norm_sq(eps_p); // = ||eps_p||^2

                            g = dot(J, eps_p); // = J' . eps_p

                            stop = std::fabs(g) < epsilon || eps_p_norm2 < epsilon * epsilon;

                            auto t = 2 * rho - 1;
                            mu *= std::max(one_three, 1 - t * t * t);
                            nu = 2;
                            break;
                        }
                        else {
                            mu *= nu;
                            nu *= 2;
                        }
                    }
                }
            }

            p0 = p;
        }
    };

    template<std::size_t N>
    struct MinimizeLM<2, N>
    {
        static_assert(N >= 2);

        // minimizing: ||x-f(p)||^2
        template<
            typename ObjectiveFunction = void(const VecMatAlg::tvec<double, 2>& p, VecMatAlg::tvec<double, N>& f),
            typename JacobianFunction  = void(const VecMatAlg::tvec<double, 2>& p, VecMatAlg::tmat<double, N>& J),
            typename VecF              = VecMatAlg::tvec<double, N>,
            typename VecP              = VecMatAlg::tvec<double, 2>
        >
        static void solve(ObjectiveFunction getF, JacobianFunction getJ, const VecF& x, VecP& p0, int max_iter = 100)
        {
            using namespace VecMatAlg;
            using type = std::remove_cv_t<std::remove_reference_t<decltype(p0[0])>>;
            using matJ = tmat<type, 2, N>; //? fix me: Jacobian should be N*2
            using vecP = tvec<type, 2>;
            using vecF = tvec<type, N>;
            using matH = tmat<type, 2>;

            constexpr auto epsilon   = std::numeric_limits<type>::epsilon();
            constexpr auto tau       = type{ 1E-3 };
            constexpr auto one_three = type{ 1 } / type{ 3 };

            auto& vec_p0 = vec_view<type, 2>(&p0[0]);
            auto& vec_x  = vec_view<type, N>(& x[0]);

            //
            type nu = 2;
            vecP p  = vec_p0;

            vecF f; // N*1
            getF(p, f);

            // J = {{df1/dp1,df1/dp2,...,df1/dpM},
            //      {df2/dp1,df2/dp2,...,df2/dpM},
            //      ...
            //      {dfN/dp1,dfN/dp2,...,dfN/dpM}}

            matJ J; // 2*N
            getJ(p, J);

            //auto JT = transpose(J); // 2*N
            matH A;// = dot(JT, J); // 2*2
            A[0][0] = dot(J[0], J[0]);
            A[1][1] = dot(J[1], J[1]);
            A[0][1] = A[1][0] = dot(J[0], J[1]);

            auto eps_p       = vec_x - f; // N*1
            auto eps_p_norm2 = norm_sq(eps_p);

            //auto g = dot(JT, eps_p); // 2*1
            vecP g = { dot(J[0],eps_p),dot(J[1],eps_p) };

            auto stop = std::fabs(g[0]) < epsilon && std::fabs(g[1]) < epsilon;
            auto mu   = tau * std::max(A[0][0], A[1][1]);

            int it = 0, n = 0;
            for (; !stop && it < max_iter; ++it) {
                for (; !stop;) {
                    ++n;

                    // solve: (A + mu I) dp = g
                    A[0][0] += mu; A[1][1] += mu;
                    inverse_symm(A);
                    auto dp = dot(A, g);

                    // stop = ||dp|| <= eps2 * ||p||
                    stop = norm_sq(dp) <= epsilon * epsilon * norm_sq(p);

                    if (!stop) {
                        // p_new = p + dp
                        auto p_new = p + dp;

                        // rho = (||eps_p||^2 - ||x-f(p)||^2)/(dp' . (mu * dp + g))
                        getF(p_new, f);
                        auto rho = (eps_p_norm2 - distance_sq(vec_x, f)) / (dot(dp, mu * dp + g));

                        if (rho > 0) {
                            p = p_new;
                            getJ(p, J);

                            //A = dot(J, J); // = J' . J
                            A[0][0] = dot(J[0], J[0]);
                            A[1][1] = dot(J[1], J[1]);
                            A[0][1] = A[1][0] = dot(J[0], J[1]);

                            eps_p       = vec_x - f;
                            eps_p_norm2 = norm_sq(eps_p); // = ||eps_p||^2

                            //g = dot(J, eps_p); // = J' . eps_p
                            g[0] = dot(J[0], eps_p);
                            g[1] = dot(J[1], eps_p);

                            stop = norm_infinity(g) < epsilon || eps_p_norm2 < epsilon * epsilon;

                            auto t = 2 * rho - 1;
                            mu *= std::max(one_three, 1 - t * t * t);
                            nu = 2;
                            break;
                        }
                        else {
                            mu *= nu;
                            nu *= 2;
                        }
                    }
                }
            }

            vec_p0 = p;
        }
    };

    template<std::size_t N>
    struct MinimizeLM<3, N>
    {
        static_assert(N >= 3);

        // minimizing: ||x-f(p)||
        template<
            typename ObjectiveFunction = void(const VecMatAlg::tvec<double, 3>& p, VecMatAlg::tvec<double, N>& f),
            typename JacobianFunction  = void(const VecMatAlg::tvec<double, 3>& p, VecMatAlg::tmat<double, N>& J),
            typename VecF              = VecMatAlg::tvec<double, N>,
            typename VecP              = VecMatAlg::tvec<double, 3>
        >
        static void solve(ObjectiveFunction getF, JacobianFunction getJ, const VecF& x, VecP& p0, int max_iter = 100)
        {
            using namespace VecMatAlg;
            using type = std::remove_cv_t<std::remove_reference_t<decltype(p0[0])>>;
            using matJ = tmat<type, 3, N>; //? fix me: Jacobian should be N*3
            using vecP = tvec<type, 3>;
            using vecF = tvec<type, N>;
            using matH = tmat<type, 3>;
            constexpr auto epsilon = std::numeric_limits<type>::epsilon();
            constexpr type tau = 1E-3;
            constexpr type one_three = type{ 1 } / type{ 3 };

            auto& vec_p0 = vec_view<type, 3>(&p0[0]);
            auto& vec_x  = vec_view<type, N>(& x[0]);

            //
            type nu = 2;
            vecP p  = vec_p0;

            matJ J; // 3*N
            vecF f; // N*1
            getF(p, f);
            getJ(p, J);

            //auto JT = transpose(J); // 3*N
            matH A; // = dot(JT, J); // 3*3
            A[0][0] = dot(J[0], J[0]); //? fix me: 
            A[1][1] = dot(J[1], J[1]);
            A[2][2] = dot(J[2], J[2]);
            A[0][1] = A[1][0] = dot(J[0], J[1]);
            A[0][2] = A[2][0] = dot(J[0], J[2]);
            A[1][2] = A[2][1] = dot(J[1], J[2]);

            auto eps_p = vec_x - f; // N*1
            auto eps_p_norm2 = norm_sq(eps_p);
            vecP g = { dot(J[0],eps_p),dot(J[1],eps_p) ,dot(J[2],eps_p) };
            auto stop = std::fabs(g[0]) < epsilon && std::fabs(g[1]) < epsilon && std::fabs(g[2]) < epsilon;
            auto mu   = tau * std::max(A[0][0], std::max(A[1][1], A[2][2]));

            int it = 0, n = 0;
            for (; !stop && it < max_iter; ++it) {
                for (; !stop;) {
                    ++n;

                    // solve: (A + mu I) dp = g
                    A[0][0] += mu; A[1][1] += mu; A[2][2] += mu;
                    inverse_symm(A);
                    auto dp = dot(A, g);

                    // stop = ||dp|| <= eps2 * ||p||
                    stop = norm_sq(dp) <= epsilon * epsilon * norm_sq(p);

                    if (!stop) {
                        // p_new = p + dp
                        auto p_new = p + dp;

                        // rho = (||eps_p||^2 - ||x-f||^2)/(dp' . (mu * dp + g))
                        getF(p_new, f);
                        auto rho = (eps_p_norm2 - distance_sq(vec_x, f)) / (dot(dp, mu * dp + g));
                        if (rho > 0) {
                            p = p_new;
                            getJ(p, J);

                            //A = dot(J, J); // = J' . J
                            A[0][0] = dot(J[0], J[0]);
                            A[1][1] = dot(J[1], J[1]);
                            A[2][2] = dot(J[2], J[2]);
                            A[0][1] = A[1][0] = dot(J[0], J[1]);
                            A[0][2] = A[2][0] = dot(J[0], J[2]);
                            A[1][2] = A[2][1] = dot(J[1], J[2]);

                            eps_p       = vec_x - f;
                            eps_p_norm2 = norm_sq(eps_p); // = ||eps_p||^2

                            //g = dot(J, eps_p); // = J' . eps_p
                            g[0] = dot(J[0], eps_p);
                            g[1] = dot(J[1], eps_p);
                            g[2] = dot(J[2], eps_p);

                            stop = norm_infinity(g) < epsilon || eps_p_norm2 < epsilon * epsilon;
                            auto t = 2 * rho - 1;
                            mu *= std::max(one_three, 1 - t * t * t);
                            nu = 2;
                            break;
                        }
                        else {
                            mu *= nu;
                            nu *= 2;
                        }
                    }
                }
            }

            vec_p0 = p;
        }
    };
}
