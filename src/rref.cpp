#include "rref.h"
#include "constraint.h"
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>


void rref(const int nrows,
          const int ncols,
          double **mat,
          int &nrank,
          const double tolerance)
{
    // Return the reduced row echelon form (rref) of matrix mat.
    // In addition, rank of the matrix is estimated.

    int jcol;
    double tmp;

    nrank = 0;

    auto icol = 0;

    for (auto irow = 0; irow < nrows; ++irow) {

        auto pivot = irow;

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

        for (auto jrow = 0; jrow < nrows; ++jrow) {
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


void rref_sparse(const int ncols,
                 ConstraintSparseForm &sp_constraint,
                 const double tolerance)
{
    // This function is somewhat sensitive to the numerical accuracy.
    // The loss of numerical digits can lead to instability.
    // Column ordering may improve the stability, but I'm not sure.
    // Smaller tolerance is preferable.

    const int nrows = sp_constraint.size();
    int jrow;
    double scaling_factor;
    double division_factor;

    // This parameter controls the stability and performance.
    // Smaller value is more stable but little more costly.
    //    double zero_criterion = tolerance * 1.0e-3;
    const auto zero_criterion = eps15;

    auto nrank = 0;
    auto icol = 0;

    std::set<unsigned int>::iterator it_found;
    std::map<unsigned int, double>::iterator it_elem;

    for (auto irow = 0; irow < nrows; ++irow) {

        auto pivot = irow;

        while (true) {
            it_elem = sp_constraint[pivot].find(icol);
            if (it_elem != sp_constraint[pivot].end()) {
                if (std::abs(it_elem->second) >= tolerance) {
                    break;
                }
            }

            ++pivot;
            if (pivot == nrows) {
                pivot = irow;
                ++icol;

                if (icol == ncols) break;
            }
        }

        if (icol == ncols) break;

        if (std::abs(it_elem->second) >= tolerance) ++nrank;

        if (pivot != irow) {
            std::iter_swap(sp_constraint.begin() + irow,
                           sp_constraint.begin() + pivot);
        }

        division_factor = 1.0 / it_elem->second;
        for (auto &it : sp_constraint[irow]) {
            it.second *= division_factor;
        }

        for (jrow = 0; jrow < nrows; ++jrow) {
            if (jrow == irow) continue;

            it_elem = sp_constraint[jrow].find(icol);
            if (it_elem == sp_constraint[jrow].end()) continue;
            scaling_factor = it_elem->second;

            // Subtract irow elements from jrow
            for (const auto &it_now : sp_constraint[irow]) {
                if (it_now.first < icol) {
                    continue;
                }
                auto it_other = sp_constraint[jrow].find(it_now.first);
                if (it_other != sp_constraint[jrow].end()) {
                    it_other->second -= scaling_factor * it_now.second;
                    // Delete zero elements and remove from map. 
                    // A smaller threshould is used for better stability.
                    if (std::abs(it_other->second) < zero_criterion) {
                        sp_constraint[jrow].erase(it_other);
                    }
                } else {
                    sp_constraint[jrow][it_now.first] = -scaling_factor * it_now.second;
                }
            }
            // Make sure to erase the icol element from the target row if it exists.
            // When the original pivot element is large, the element after subtraction can sometimes be 
            // larger than the tolerance value because of the loss of significant digis.
            auto it_other = sp_constraint[jrow].find(icol);
            if (it_other != sp_constraint[jrow].end()) {
                sp_constraint[jrow].erase(it_other);
            }
        }
    }

    // Erase all elements smaller than the tolerance value
    for (jrow = 0; jrow < nrows; ++jrow) {
        auto it_other = sp_constraint[jrow].begin();
        while (it_other != sp_constraint[jrow].end()) {
            if (std::abs(it_other->second) <= tolerance) {
                sp_constraint[jrow].erase(it_other++);
            } else {
                ++it_other;
            }
        }
    }

    // Remove emptry entries from the sp_constraint vector
    sp_constraint.erase(std::remove_if(sp_constraint.begin(),
                                       sp_constraint.end(),
                                       [](const std::map<unsigned int, double> &obj) { return obj.empty(); }),
                        sp_constraint.end());
    sp_constraint.shrink_to_fit();
}
