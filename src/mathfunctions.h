/*
 mathfunctions.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <iostream>
#include <cstdlib>
#include <vector>

template <typename T>
inline void matmul3(T ret[3][3], const T amat[3][3], const T bmat[3][3]) {
    int i, j, k;

    T ret_tmp[3][3];

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ret_tmp[i][j] = 0.0;
            for (k = 0; k < 3; ++k) ret_tmp[i][j] += amat[i][k] * bmat[k][j]; 	        
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ret[i][j] = ret_tmp[i][j];
        }
    }
}

inline void transpose3(double ret[3][3], const double mat[3][3]) 
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ret[i][j] = mat[j][i];
        }
    }
}

inline void rotvec(double vec_out[3], double vec_in[3], double mat[3][3], char mode = 'N')
{
    // Perform matrix x vector multiplication. 
    //
    // vec_out = mat      * vec_in   (mode = 'N')
    //          (mat)^{t} * vec_in   (mode = 'T')
    //

    unsigned int i;
    double vec_tmp[3];

    for (i = 0; i < 3; ++i){
        vec_tmp[i] = vec_in[i];
    }

    if (mode == 'N') {
        for (i = 0; i < 3; ++i){
            vec_out[i] = mat[i][0] * vec_tmp[0] + mat[i][1] * vec_tmp[1] + mat[i][2] * vec_tmp[2];
        }
    } else if (mode == 'T'){
        for (i = 0; i < 3; ++i){
            vec_out[i] = mat[0][i] * vec_tmp[0] + mat[1][i] * vec_tmp[1] + mat[2][i] * vec_tmp[2];
        }
    } else {
        std::cout << "Invalid mode " << mode << std::endl;
        exit(1);
    }
}

inline void rotvec(double vec_out[3], double vec_in[3], double **mat, char mode = 'N')
{
    // Perform matrix x vector multiplication. 
    //
    // vec_out = mat      * vec_in   (mode = 'N')
    //          (mat)^{t} * vec_in   (mode = 'T')
    //

    unsigned int i;
    double vec_tmp[3];

    for (i = 0; i < 3; ++i){
        vec_tmp[i] = vec_in[i];
    }

    if (mode == 'N') {
        for (i = 0; i < 3; ++i){
            vec_out[i] = mat[i][0] * vec_tmp[0] + mat[i][1] * vec_tmp[1] + mat[i][2] * vec_tmp[2];
        }
    } else if (mode == 'T'){
        for (i = 0; i < 3; ++i){
            vec_out[i] = mat[0][i] * vec_tmp[0] + mat[1][i] * vec_tmp[1] + mat[2][i] * vec_tmp[2];
        }
    } else {
        std::cout << "Invalid mode " << mode << std::endl;
        exit(1);
    }
}

inline void invmat3(double invmat[3][3], double mat[3][3])
{
    unsigned int i, j;
    double det;
    double mat_tmp[3][3];

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            mat_tmp[i][j] = mat[i][j];
        }
    }

    det = mat_tmp[0][0] * mat_tmp[1][1] * mat_tmp[2][2] 
    + mat_tmp[1][0] * mat_tmp[2][1] * mat_tmp[0][2] 
    + mat_tmp[2][0] * mat_tmp[0][1] * mat_tmp[1][2]
    - mat_tmp[0][0] * mat_tmp[2][1] * mat_tmp[1][2] 
    - mat_tmp[2][0] * mat_tmp[1][1] * mat_tmp[0][2]
    - mat_tmp[1][0] * mat_tmp[0][1] * mat_tmp[2][2];

    if(std::abs(det) < 1.0e-12) {
        std::cout << "invmat3: Given matrix is singular" << std::endl;
        exit(1);
    }

    double factor = 1.0 / det;

    invmat[0][0] = (mat_tmp[1][1] * mat_tmp[2][2] - mat_tmp[1][2] * mat_tmp[2][1]) * factor;
    invmat[0][1] = (mat_tmp[0][2] * mat_tmp[2][1] - mat_tmp[0][1] * mat_tmp[2][2]) * factor;
    invmat[0][2] = (mat_tmp[0][1] * mat_tmp[1][2] - mat_tmp[0][2] * mat_tmp[1][1]) * factor;

    invmat[1][0] = (mat_tmp[1][2] * mat_tmp[2][0] - mat_tmp[1][0] * mat_tmp[2][2]) * factor;
    invmat[1][1] = (mat_tmp[0][0] * mat_tmp[2][2] - mat_tmp[0][2] * mat_tmp[2][0]) * factor;
    invmat[1][2] = (mat_tmp[0][2] * mat_tmp[1][0] - mat_tmp[0][0] * mat_tmp[1][2]) * factor;

    invmat[2][0] = (mat_tmp[1][0] * mat_tmp[2][1] - mat_tmp[1][1] * mat_tmp[2][0]) * factor;
    invmat[2][1] = (mat_tmp[0][1] * mat_tmp[2][0] - mat_tmp[0][0] * mat_tmp[2][1]) * factor;
    invmat[2][2] = (mat_tmp[0][0] * mat_tmp[1][1] - mat_tmp[0][1] * mat_tmp[1][0]) * factor;
}

inline void invmat3_i(int invmat[3][3], int mat[3][3])
{
    int det;

    det = mat[0][0] * mat[1][1] * mat[2][2] 
    + mat[1][0] * mat[2][1] * mat[0][2] 
    + mat[2][0] * mat[0][1] * mat[1][2]
    - mat[0][0] * mat[2][1] * mat[1][2] 
    - mat[2][0] * mat[1][1] * mat[0][2]
    - mat[1][0] * mat[0][1] * mat[2][2];

    if(std::abs(det) == 0) {
        std::cout << "invmat3_i: Given matrix is singular" << std::endl;
        exit(1);
    }

    invmat[0][0] = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) / det;
    invmat[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) / det;
    invmat[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) / det;

    invmat[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) / det;
    invmat[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) / det;
    invmat[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) / det;

    invmat[2][0] = (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) / det;
    invmat[2][1] = (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) / det;
    invmat[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) / det;

}

inline int nint(double x)
{
    return int(x + 0.5 - (x < 0.0));
}

template <typename T>
void insort(int n, T *arr)
{
    int i, j;
    T tmp;

    for (i = 1; i < n; ++i) {
        tmp = arr[i];
        for (j = i - 1; j >= 0 && arr[j] > tmp; --j) {
            arr[j + 1] = arr[j];
        }
        arr[j + 1] = tmp;
    }
}

inline void sort_tail(const int n, int *arr)
{
    int i, m;

    m = n - 1;
    int *ind_tmp;

    ind_tmp = new int[m];

    for (i = 0; i < m; ++i) {
        ind_tmp[i] = arr[i + 1];
    }

    insort(m, ind_tmp);

    for (i = 0; i < m; ++i) {
        arr[i + 1] = ind_tmp[i];
    }
    delete [] ind_tmp;
}


void rref(const int nrows,
    const int ncols,
    double **mat,
    int &nrank,
    const double tolerance)
{
    // Return the reduced row echelon form (rref) of matrix mat.
    // In addition, rank of the matrix is estimated.

    int irow, icol, jrow, jcol;
    int pivot;
    double tmp;

    nrank = 0;

    icol = 0;

    for (irow = 0; irow < nrows; ++irow) {

        pivot = irow;

        while (std::abs(mat[pivot][icol]) < tolerance) {
            ++pivot;

            if (pivot == nrows) {
                pivot = irow;
                ++icol;

                if (icol == ncols) break;
            }
        }

        if (icol == ncols) break;

        if (std::abs(mat[pivot][icol]) > tolerance) ++nrank;

        if (pivot != irow) {
            //#pragma omp parallel for private(tmp)
            for (jcol = icol; jcol < ncols; ++jcol) {
                tmp = mat[pivot][jcol];
                mat[pivot][jcol] = mat[irow][jcol];
                mat[irow][jcol] = tmp;
            }
        }

        tmp = mat[irow][icol];
        tmp = 1.0 / tmp;
        //#pragma omp parallel for
        for (jcol = icol; jcol < ncols; ++jcol) {
            mat[irow][jcol] *= tmp;
        }

        for (jrow = 0; jrow < nrows; ++jrow) {
            if (jrow == irow) continue;

            tmp = mat[jrow][icol];
            //#pragma omp parallel for
            for (jcol = icol; jcol < ncols; ++jcol) {
                mat[jrow][jcol] -= tmp * mat[irow][jcol];
            }
        }
    }
}


void rref(std::vector<std::vector<double>> &mat,
    const double tolerance)
{
    // Return the reduced row echelon form (rref) of matrix mat.
    // In addition, rank of the matrix is estimated.

    int irow, icol, jrow, jcol;
    int pivot;
    double tmp;

    int nrank = 0;

    icol = 0;

    int nrows = mat.size();
    int ncols = mat[0].size();

    for (irow = 0; irow < nrows; ++irow) {

        pivot = irow;

        while (std::abs(mat[pivot][icol]) < tolerance) {
            ++pivot;

            if (pivot == nrows) {
                pivot = irow;
                ++icol;

                if (icol == ncols) break;
            }
        }

        if (icol == ncols) break;

        if (std::abs(mat[pivot][icol]) > tolerance) ++nrank;

        if (pivot != irow) {
            for (jcol = icol; jcol < ncols; ++jcol) {
                tmp = mat[pivot][jcol];
                mat[pivot][jcol] = mat[irow][jcol];
                mat[irow][jcol] = tmp;
            }
        }

        tmp = mat[irow][icol];
        tmp = 1.0 / tmp;
        for (jcol = icol; jcol < ncols; ++jcol) {
            mat[irow][jcol] *= tmp;
        }

        for (jrow = 0; jrow < nrows; ++jrow) {
            if (jrow == irow) continue;

            tmp = mat[jrow][icol];
            for (jcol = icol; jcol < ncols; ++jcol) {
                mat[jrow][jcol] -= tmp * mat[irow][jcol];
            }
        }
    }

    mat.erase(mat.begin() + nrank, mat.end());
    mat.shrink_to_fit();
}

