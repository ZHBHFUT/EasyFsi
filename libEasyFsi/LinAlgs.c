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

//!-------------------------------------------------------------
//! @file       LinAlgs.c
//!             The implement of linear algebraic functions.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------
#include "LinAlgs.h"

#include <math.h>
#include <stdlib.h>

#include "Assert.hpp"

void mat_inverse(int rank, double* m, int* buffer, int* singular)
{
    int i, j, k, l, u, v, *is, *js, alloc = 0;
    double d, p;

    *singular = 0;

    if (rank == 0)return;

    if (rank == 1) {
        m[0] = 1.0 / m[0];
        *singular = m[0] == 0;
        return;
    }

    // 求出 m[rank][rank] 之逆矩阵，并将所求逆矩阵置入 m[rank][rank]
    // 全选主元 Gauss-Jordan 消去法

    if (!buffer) {
        buffer = (int*)malloc(sizeof(int) * rank * 2);
        alloc = 1;
    }

    is = buffer;        // size=rank
    js = buffer + rank; // size=rank

    if (!is || !js) { ASSERT(0); abort(); }

    for (k = 0; k < rank; k++) {
        d = 0;

        for (i = k; i < rank; i++) {
            for (j = k; j < rank; j++) {
                l = i * rank + j;
                p = fabs(m[l]);

                if (p > d) {
                    d = p;
                    is[k] = i;
                    js[k] = j;
                }
            }
        }

        if (d == 0) {
            ASSERT(0);
            *singular = 1;
            if (alloc)free(buffer);
            return;
        }

        if (is[k] != k) {
            for (j = 0; j < rank; j++) {
                u = k * rank + j;
                v = is[k] * rank + j;
                p = m[u]; m[u] = m[v]; m[v] = p;
            }
        }

        if (js[k] != k) {
            for (i = 0; i < rank; i++) {
                u = i * rank + k;
                v = i * rank + js[k];

                p = m[u]; m[u] = m[v]; m[v] = p;
            }
        }

        l = k * rank + k;
        m[l] = 1 / m[l];
        for (j = 0; j < rank; j++) {
            if (j != k) {
                u = k * rank + j;
                m[u] = m[u] * m[l];
            }
        }

        for (i = 0; i < rank; i++) {
            if (i != k) {
                for (j = 0; j < rank; j++) {
                    if (j != k) {
                        u = i * rank + j;
                        m[u] = m[u] - m[i * rank + k] * m[k * rank + j];
                    }
                }
            }
        }
        for (i = 0; i < rank; i++) {
            if (i != k) {
                u = i * rank + k;
                m[u] = -m[u] * m[l];
            }
        }
    }

    for (k = rank - 1; k >= 0; k--) {
        if (js[k] != k) {
            for (j = 0; j < rank; j++) {
                u = k * rank + j;
                v = js[k] * rank + j;

                p = m[u]; m[u] = m[v]; m[v] = p;
            }
        }
        if (is[k] != k) {
            for (i = 0; i < rank; i++) {
                u = i * rank + k;
                v = i * rank + is[k];
                p = m[u]; m[u] = m[v]; m[v] = p;
            }
        }
    }

    if (alloc)free(buffer);
}
void mat_transpose(int rank, double* m)
{
    int i, j;
    double t;
    for (i = 0; i < rank; ++i) {
        for (j = i + 1; j < rank; ++j) {
            //std::swap(m[i * rank + j], m[j * rank + i]);
            t = m[i * rank + j];
            m[i * rank + j] = m[j * rank + i];
            m[j * rank + i] = t;
        }
    }
}
void mat_transpose_to(int m, int n, const double* A, double* AT)
{
    ASSERT(A != AT);
    int i, j;
    for (i = 0; i < n; ++i) {
        for (j = i + 1; j < m; ++j, ++AT) {
            *AT = A[j * n + i];
        }
    }
}
void mat_apply_vec(int nrow, int ncol, const double* A, const double* x, double* y)
{
    ASSERT(x != y);
    int i, j;
    const double* xi;
    for (i = 0; i < nrow; ++i, ++y) {
        *y = 0;
        xi = x;
        for (j = 0; j < ncol; ++j, ++A, ++xi) {
            (*y) += (*A) * (*xi);
        }
    }
}
void mat_transpose_apply_vec(int nrow, int ncol, const double* A, const double* x, double* y)
{
    ASSERT(x != y);
    int i, j;
    const double* r;
    for (i = 0; i < ncol; ++i, ++y) {
        //y[i] = dot(A(:,i),x)
        *y = 0;
        r = x;
        for (j = 0; j < nrow; ++j, ++r) {
            (*y) += A[j * ncol + i] * (*r);
        }
    }
}
void mat_apply_mat(int m, int n, int l, const double* A, const double* B, double* C)
{
    ASSERT(A != B && B != C);
    int i, j, k;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < l; ++j, ++C) {
            // C(i,j) = A(i,:).B(:,j)
            *C = 0;
            for (k = 0; k < n; ++k) {
                *C += A[k] * B[k * l + j];
            }
        }
        A += n;
    }
}
void mat_transpose_apply_mat(int m, int n, int l, const double* A, const double* x, double* y)
{
    ASSERT(x != y);
    int i, j, k;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < l; ++j, ++y) {
            // y(i,j) = dot(A(:,i),x(:,j)
            *y = 0;
            for (k = 0; k < m; ++k) {
                *y += A[k * n + i] * x[k * l + j];
            }
        }
    }
}
void mat_apply_add_vec(int nrow, int ncol, const double* A, const double* x, double* y)
{
    ASSERT(x != y);
    int i, j;
    const double* xi;
    double t;
    for (i = 0; i < nrow; ++i, ++y) {
        xi = x;
        t = 0;
        for (j = 0; j < ncol; ++j, ++A, ++xi) {
            t += (*A) * (*xi);
        }
        *y += t;
    }
}
void mat_apply_add_mat(int m, int n, int l, const double* A, const double* B, double* C)
{
    ASSERT(A != B && B != C);
    int i, j, k;
    double t;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < l; ++j, ++C) {
            t = 0;
            for (k = 0; k < n; ++k) {
                t += A[k] * B[k * l + j];
            }
            *C += t;
        }
        A += n;
    }
}
void compute_xps_ibuffer_size(int np_src, int ndim, int* ibuffer_size)
{
    int rank = np_src + 1 + ndim;
    *ibuffer_size = rank * 2;
}
void compute_xps_dbuffer_size(int np_src, int ndim, int* dbuffer_size)
{
    int rank = np_src + 1 + ndim;
    *dbuffer_size = rank * rank + rank;
}
void compute_xps_interp_matrix(int np_src, int np_des, int ndim, int coord_stride, const double* coords_src, const double* coords_des, double* mat_G, double* dbuffer, int* ibuffer, int* singular)
{
#define eps  1.0E-40
    const int rank = np_src + 1 + ndim;// size of matrix TS

    double dx, rij2, m, *TS, *TSi, *TFi;
    const double* p, * q;
    int i, j, k, alloc = 0;
    
    if (!dbuffer) {
        dbuffer = (double*)malloc(sizeof(double) * (rank * rank + rank));
        alloc = 1;
    }

    TS  = dbuffer;
    TFi = dbuffer + rank * rank;

    if (!TS || !TFi || !ibuffer) { ASSERT(0); }

    //--------------------------------------------
    // Step 1: generate TS matrix
    //--------------------------------------------
    // 
    //   r00 r01 ... r0N  1  x0  y0  z0
    //   r10 r11 ... r1N  1  x1  y1  z1
    //   ...
    //   rN0 rN1 ... rNN  1  xN  yN  zN
    //     1   1 ...   1  0   0   0   0
    //    x0  x1 ...  xN  0   0   0   0
    //    y0  y1 ...  yN  0   0   0   0
    //    z0  z1 ...  zN  0   0   0   0
    //

    // 1) left-top block: r^2 log(r^2)
    TSi = TS;
    for (i = 0, p = coords_src; i < np_src; ++i) {
        for (j = 0, q = coords_src; j < np_src; ++j) {
            // distance
            rij2 = 0;
            for (k = 0; k < ndim; ++k) {
                dx = p[k] - q[k];
                rij2 += dx * dx;
            }
            TSi[j] = rij2 * log(rij2 + eps);// prevent float overflow

            q += coord_stride;
        }
        TSi += rank;
        p += coord_stride;
    }
    // 2) right top block: 1, x, y, z
    for (i = 0, p = coords_src; i < np_src; ++i) {
        TS[i * rank + np_src + 0] = 1;
        for (j = 0; j < ndim; ++j)
            TS[i * rank + np_src + j + 1] = p[j];
        p += coord_stride;
    }
    // 3) left bottom block: 
    //   1   1   ...
    //   x0  x1  ...
    //   y0  y1  ...
    //   .   .   ...
    for (j = 0; j < np_src; ++j) {
        TS[(np_src + 0) * rank + j] = 1;
        for (k = 0; k < ndim; ++k)
            TS[(np_src + k + 1) * rank + j] = coords_src[coord_stride * j + k];
    }
    // 4) right bottom block: zeros
    for (i = np_src; i < rank; ++i) {
        for (j = np_src; j < rank; j++)
            TS[i * rank + j] = 0;
    }

    //--------------------------------------------
    // Step 2: inverse TS matrix
    //--------------------------------------------

    mat_inverse(rank, TS, ibuffer, singular);
    if (*singular) {
        if (alloc)free(dbuffer);
        return;
    }

    //-------------------------------------------------------------------------
    // Step 3: G = Tf.Ts^-1
    //-------------------------------------------------------------------------
    // NOTE: only output left-top block of G

    // TF = 
    //   r00  r01  ... r0N  1  x0  y0  z0
    //   r10  r11  ... r1N  1  x1  y1  z1
    //   ...
    //   rM0  rM1  ... rMN  1  xM  yM  zM
    //
    
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
            TFi[j] = rij2 * log(rij2 + eps);

            q += coord_stride;
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

        p += coord_stride;
    }

    if (alloc)free(dbuffer);
#undef eps
}
