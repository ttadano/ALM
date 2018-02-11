#include "rref.h"
#include "memory.h"
#include <vector>

void remove_redundant_rows(const int n,
                           std::vector<std::vector<double>> &constraint_mat,
                           const double tolerance)
{
    int i, j;

    int nparam = n;
    int nconst = constraint_mat.size();
    double **mat_tmp;
    std::vector<double> arr_tmp;

    int nrank;

    if (nconst > 0) {

        allocate(mat_tmp, nconst, nparam);

        for (i = 0; i < nconst; ++i) {
            for (j = 0; j < nparam; ++j) {
                mat_tmp[i][j] = constraint_mat[i][j];
            }
        }

        rref(nconst, nparam, mat_tmp, nrank, tolerance);
        arr_tmp.resize(nparam);
        constraint_mat.clear();

        for (i = 0; i < nrank; ++i) {
            for (j = 0; j < i; ++j) arr_tmp[j] = 0.0;

            for (j = i; j < nparam; ++j) {
                if (std::abs(mat_tmp[i][j]) < tolerance) {
                    arr_tmp[j] = 0.0;
                } else {
                    arr_tmp[j] = mat_tmp[i][j];
                }
            }
            constraint_mat.push_back(arr_tmp);
        }

        deallocate(mat_tmp);

    }
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
