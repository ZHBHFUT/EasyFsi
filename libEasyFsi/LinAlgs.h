#pragma once

#ifdef __cplusplus
extern "C" {
#endif

    //-------------------------------------------------
    // ���Դ������㺯��
    //-------------------------------------------------

    //! @brief ������� m �������
    //! @param rank ����������
    //! @param m    ����Ԫ���׵�ַ
    //! @param buffer ������������Ϊ 2*rank�����Ϊnull�����ڲ���̬���䡣
    //! @param singular ���������0=��������죬1=��������
    void mat_inverse(int rank, double* m, int* buffer, int* singular);

    //! @brief �Է������ת�ò���
    //! @param rank ����������
    //! @param m ����Ԫ���׵�ַ
    void mat_transpose(int rank, double* m);

    //! @brief ����ת��
    //! @param m ��������
    //! @param n ��������
    //! @param A ���վ����Ԫ���׵�ַ
    //! @param AT ��������Ԫ���׵�ַ
    void mat_transpose_to(int m, int n, const double* A, double* AT);

    //! @brief ������������˻��� y = A * x
    //! @param nrow ��������
    //! @param ncol ��������
    //! @param A ����Ԫ���׵�ַ
    //! @param x ����xԪ���׵�ַ
    //! @param y �������y��Ԫ���׵�ַ
    //! @note ����\x��Ԫ�ظ���Ӧ��Ϊ\ncol��\y��Ԫ�ظ���Ӧ��Ϊ\nrow.
    void mat_apply_vec(int nrow, int ncol, const double* A, const double* x, double* y);

    //! @brief �������ת���������ĳ˻��� y = A^T * x
    //! @param nrow ��������
    //! @param ncol ��������
    //! @param A ����Ԫ���׵�ַ
    //! @param x ����xԪ���׵�ַ
    //! @param y �������y��Ԫ���׵�ַ
    //! @note ����\x��Ԫ�ظ���Ӧ��Ϊ\nrow��\y��Ԫ�ظ���Ӧ��Ϊ\ncol.
    void mat_transpose_apply_vec(int nrow, int ncol, const double* A, const double* x, double* y);

    //! @brief ������������ĳ˻��� C = A * B
    //! @param m ����A������
    //! @param n ����A������
    //! @param l ����B������
    //! @param A ����AԪ���׵�ַ��Ԫ�ظ���Ϊ m*n
    //! @param B ����BԪ���׵�ַ��Ԫ�ظ���Ϊ n*l
    //! @param C �������C��Ԫ���׵�ַ��Ԫ�ظ���Ϊ m*l
    void mat_apply_mat(int m, int n, int l, const double* A, const double* B, double* C);

    //! @brief �������ת���������ĳ˻��� y = A^T * x
    //! @param m ����A������
    //! @param n ����A������
    //! @param l ����B������
    //! @param A ����Ԫ���׵�ַ���ߴ�Ϊ m*n
    //! @param x ����xԪ���׵�ַ���ߴ�Ϊ m*l
    //! @param y �������y��Ԫ���׵�ַ���ߴ�Ϊ n*l
    void mat_transpose_apply_mat(int m, int n, int l, const double* A, const double* x, double* y);

    //! @brief ������������˻��� y = A * x + y
    //! @param nrow ��������
    //! @param ncol ��������
    //! @param A ����Ԫ���׵�ַ���ߴ�Ϊ m*n
    //! @param x ����xԪ���׵�ַ���ߴ�Ϊ ncol
    //! @param y �������y��Ԫ���׵�ַ���ߴ�Ϊ nrow
    void mat_apply_add_vec(int nrow, int ncol, const double* A, const double* x, double* y);

    //! @brief ������������ĳ˻��� C = A * B + C
    //! @param m ����A������
    //! @param n ����A������
    //! @param l ����B������
    //! @param A ����AԪ���׵�ַ��Ԫ�ظ���Ϊ m*n
    //! @param B ����BԪ���׵�ַ��Ԫ�ظ���Ϊ n*l
    //! @param C �������C��Ԫ���׵�ַ��Ԫ�ظ���Ϊ m*l
    void mat_apply_add_mat(int m, int n, int l, const double* A, const double* B, double* C);

    //-------------------------------------------------
    // ��ֵ������㺯��
    //-------------------------------------------------
    
    //! @brief ����������ֵ���������������������С
    //! @param np_src ��֪�����
    //! @param ndim   �������������1=spline, 2=ips, 3=tps
    //! @param ibuffer_size ��������������
    void compute_xps_ibuffer_size(int np_src, int ndim, int* ibuffer_size);

    //! @brief ����������ֵ��������ĸ�������������С
    //! @param np_src ��֪�����
    //! @param ndim   �������������1=spline, 2=ips, 3=tps
    //! @param dbuffer_size ����������������
    void compute_xps_dbuffer_size(int np_src, int ndim, int* dbuffer_size);

    //! @brief ����������ֵ����
    //! @param [in]     np_src        ��֪�����
    //! @param [in]     np_des        δ֪�㣨����ֵ�㣩����
    //! @param [in]     ndim          �������������1=spline, 2=ips, 3=tps
    //! @param [in]     coord_stride  ��������Ĳ�����ͨ���� \ndim ���
    //! @param [in]     coords_src    ��֪���������飬size=np_src*coord_stride���洢��ʽΪ: [x0,y0,...,x1,y1,...,]
    //! @param [in]     coords_des    δ֪���������飬size=np_des*coord_stride���洢��ʽΪ: [x0,y0,...,x1,y1,...,]
    //! @param [out]    mat_G         ���ز�ֵ���󣬳ߴ�Ϊ np_des * np_src
    //! @param [in,out] dbuffer       �������������������� \compute_xps_dbuffer_size ��ȡ�����Ϊnull�����ڲ���̬����
    //! @param [in,out] ibuffer       ������������������ \compute_xps_ibuffer_size ��ȡ�����Ϊnull�����ڲ���̬����
    //! @param [out]singular          �����Ƿ����죬0=��other=��
    void compute_xps_interp_matrix(int np_src, int np_des, int ndim, int coord_stride, const double* coords_src, const double* coords_des, double* mat_G, double* dbuffer, int* ibuffer, int* singular);
    
#ifdef __cplusplus
}
#endif
