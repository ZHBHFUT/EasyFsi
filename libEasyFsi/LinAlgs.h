#pragma once

#ifdef __cplusplus
extern "C" {
#endif

    //-------------------------------------------------
    // 线性代数运算函数
    //-------------------------------------------------

    //! @brief 计算矩阵 m 的逆矩阵
    //! @param rank 矩阵行列数
    //! @param m    矩阵元素首地址
    //! @param buffer 缓冲区，长度为 2*rank，如果为null则在内部动态分配。
    //! @param singular 输出参数，0=矩阵非奇异，1=矩阵奇异
    void mat_inverse(int rank, double* m, int* buffer, int* singular);

    //! @brief 对方阵进行转置操作
    //! @param rank 矩阵行列数
    //! @param m 矩阵元素首地址
    void mat_transpose(int rank, double* m);

    //! @brief 矩阵转置
    //! @param m 矩阵行数
    //! @param n 矩阵列数
    //! @param A 数日矩阵的元素首地址
    //! @param AT 输出矩阵的元素首地址
    void mat_transpose_to(int m, int n, const double* A, double* AT);

    //! @brief 计算矩阵向量乘积： y = A * x
    //! @param nrow 矩阵行数
    //! @param ncol 矩阵列数
    //! @param A 矩阵元素首地址
    //! @param x 向量x元素首地址
    //! @param y 输出向量y的元素首地址
    //! @note 向量\x的元素个数应该为\ncol，\y的元素个数应该为\nrow.
    void mat_apply_vec(int nrow, int ncol, const double* A, const double* x, double* y);

    //! @brief 计算矩阵转置与向量的乘积： y = A^T * x
    //! @param nrow 矩阵行数
    //! @param ncol 矩阵列数
    //! @param A 矩阵元素首地址
    //! @param x 向量x元素首地址
    //! @param y 输出向量y的元素首地址
    //! @note 向量\x的元素个数应该为\nrow，\y的元素个数应该为\ncol.
    void mat_transpose_apply_vec(int nrow, int ncol, const double* A, const double* x, double* y);

    //! @brief 计算矩阵与矩阵的乘积： C = A * B
    //! @param m 矩阵A的行数
    //! @param n 矩阵A的列数
    //! @param l 矩阵B的列数
    //! @param A 矩阵A元素首地址，元素个数为 m*n
    //! @param B 矩阵B元素首地址，元素个数为 n*l
    //! @param C 输出矩阵C的元素首地址，元素个数为 m*l
    void mat_apply_mat(int m, int n, int l, const double* A, const double* B, double* C);

    //! @brief 计算矩阵转置与向量的乘积： y = A^T * x
    //! @param m 矩阵A的行数
    //! @param n 矩阵A的列数
    //! @param l 矩阵B的列数
    //! @param A 矩阵元素首地址，尺寸为 m*n
    //! @param x 矩阵x元素首地址，尺寸为 m*l
    //! @param y 输出矩阵y的元素首地址，尺寸为 n*l
    void mat_transpose_apply_mat(int m, int n, int l, const double* A, const double* x, double* y);

    //! @brief 计算矩阵向量乘积： y = A * x + y
    //! @param nrow 矩阵行数
    //! @param ncol 矩阵列数
    //! @param A 矩阵元素首地址，尺寸为 m*n
    //! @param x 向量x元素首地址，尺寸为 ncol
    //! @param y 输出向量y的元素首地址，尺寸为 nrow
    void mat_apply_add_vec(int nrow, int ncol, const double* A, const double* x, double* y);

    //! @brief 计算矩阵与矩阵的乘积： C = A * B + C
    //! @param m 矩阵A的行数
    //! @param n 矩阵A的列数
    //! @param l 矩阵B的列数
    //! @param A 矩阵A元素首地址，元素个数为 m*n
    //! @param B 矩阵B元素首地址，元素个数为 n*l
    //! @param C 输出矩阵C的元素首地址，元素个数为 m*l
    void mat_apply_add_mat(int m, int n, int l, const double* A, const double* B, double* C);

    //-------------------------------------------------
    // 插值矩阵计算函数
    //-------------------------------------------------
    
    //! @brief 计算样条插值矩阵所需的整数缓冲区大小
    //! @param np_src 已知点个数
    //! @param ndim   坐标分量个数：1=spline, 2=ips, 3=tps
    //! @param ibuffer_size 整数缓冲区长度
    void compute_xps_ibuffer_size(int np_src, int ndim, int* ibuffer_size);

    //! @brief 计算样条插值矩阵所需的浮点数缓冲区大小
    //! @param np_src 已知点个数
    //! @param ndim   坐标分量个数：1=spline, 2=ips, 3=tps
    //! @param dbuffer_size 浮点数缓冲区长度
    void compute_xps_dbuffer_size(int np_src, int ndim, int* dbuffer_size);

    //! @brief 计算样条插值矩阵
    //! @param [in]     np_src        已知点个数
    //! @param [in]     np_des        未知点（待插值点）个数
    //! @param [in]     ndim          坐标分量个数：1=spline, 2=ips, 3=tps
    //! @param [in]     coord_stride  坐标数组的步长，通常与 \ndim 相等
    //! @param [in]     coords_src    已知点坐标数组，size=np_src*coord_stride，存储格式为: [x0,y0,...,x1,y1,...,]
    //! @param [in]     coords_des    未知点坐标数组，size=np_des*coord_stride，存储格式为: [x0,y0,...,x1,y1,...,]
    //! @param [out]    mat_G         返回插值矩阵，尺寸为 np_des * np_src
    //! @param [in,out] dbuffer       浮点数缓冲区，长度由 \compute_xps_dbuffer_size 获取，如果为null则在内部动态分配
    //! @param [in,out] ibuffer       整数缓冲区，长度由 \compute_xps_ibuffer_size 获取，如果为null则在内部动态分配
    //! @param [out]singular          矩阵是否奇异，0=否，other=是
    void compute_xps_interp_matrix(int np_src, int np_des, int ndim, int coord_stride, const double* coords_src, const double* coords_des, double* mat_G, double* dbuffer, int* ibuffer, int* singular);
    
#ifdef __cplusplus
}
#endif
