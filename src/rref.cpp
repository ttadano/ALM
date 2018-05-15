#include "rref.h"
#include "memory.h"
#include "constraint.h"
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

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
            for (j = 0; j < nparam; ++j) arr_tmp[j] = 0.0;
            int iloc = -1;
            for (j = 0; j < nparam; ++j) {
                if (std::abs(mat_tmp[i][j]) < tolerance) {
                    arr_tmp[j] = 0.0;
                } else {
                    arr_tmp[j] = mat_tmp[i][j];
                }

                if (std::abs(arr_tmp[j]) >= tolerance) {
                    iloc = j;
                }
            }
            if (iloc != -1) {
                constraint_mat.push_back(arr_tmp);
            }
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


void rref_sparse(const int ncols,
                 ConstraintSparseForm &sp_constraint,
                 const double tolerance)
{
    std::map<unsigned int, double> const_tmp2;
    int nconst = sp_constraint.size();
    int icol, irow, jrow;
    double scaling_factor;
    double division_factor;

    std::vector<std::vector<unsigned int>> FirstColList;
    FirstColList.resize(ncols);

    for (irow = 0; irow < nconst; ++irow) {
        FirstColList[sp_constraint[irow].begin()->first].push_back(irow);
    }

    unsigned int nsize_firstcol;

    for (icol = 0; icol < ncols; ++icol) {

        nsize_firstcol = FirstColList[icol].size();

        if (nsize_firstcol == 0) continue;

        irow = FirstColList[icol][0];
        auto const_now = sp_constraint[irow];


        division_factor = sp_constraint[irow].begin()->second;
        division_factor = 1.0 / division_factor;

        for (auto &it : sp_constraint[irow]) {
            it.second *= division_factor;
        }

        if (nsize_firstcol >= 2) {

            // Loop over rows with the same FirstCols and do reduction
            // jrow > irow
            for (auto it = FirstColList[icol].rbegin(); it != FirstColList[icol].rend() - 1; ++it) {
                jrow = (*it);

                // First non-zero element of the jrow
                scaling_factor = sp_constraint[jrow].find(icol)->second;

                // Subtract irow elements from jrow
                for (const auto &it_now : sp_constraint[irow]) {
                    auto it_other = sp_constraint[jrow].find(it_now.first);
                    if (it_other == sp_constraint[jrow].end()) {
                        sp_constraint[jrow][it_now.first] = -scaling_factor * it_now.second;
                    } else {
                        it_other->second -= scaling_factor * it_now.second;
                        //if (std::abs(it_other->second) < tolerance) {
                        //    sp_constraint[jrow].erase(it_other);
                        //}
                    }
                }
                for (auto it_other = sp_constraint[jrow].begin(); it_other != sp_constraint[jrow].end(); ++it_other) {
                    if (std::abs(it_other->second) < tolerance) {
                        sp_constraint[jrow].erase(it_other);
                    }
                }

                if (!sp_constraint[jrow].empty()) {
                    scaling_factor = sp_constraint[jrow].begin()->second;
                    if (std::abs(scaling_factor - 1.0) > tolerance) {
                        for (auto &it_other : sp_constraint[jrow]) {
                            it_other.second /= scaling_factor;
                        }
                    }
                    FirstColList[sp_constraint[jrow].begin()->first].push_back(jrow);
                }
            }
        } // close if (nsize_firstcol >= 2)

        // Subtract elements from other rows (rref)
        for (jrow = 0; jrow < nconst; ++jrow) {
            if (jrow == irow) continue;

            // Iterator of first element of the jrow
            auto it_other = sp_constraint[jrow].find(icol);
            if (it_other == sp_constraint[jrow].end()) continue;
            // First element of the jrow
            scaling_factor = it_other->second;

            // Subtract irow elements from jrow
            for (const auto &it_now : sp_constraint[irow]) {
                auto it_other = sp_constraint[jrow].find(it_now.first);
                if (it_other != sp_constraint[jrow].end()) {
                    it_other->second -= scaling_factor * it_now.second;
                } else {
                    sp_constraint[jrow][it_now.first] = -scaling_factor * it_now.second;
                }
            }
            for (auto it_other = sp_constraint[jrow].begin(); it_other != sp_constraint[jrow].end(); ++it_other) {
                if (std::abs(it_other->second) < tolerance) {
                    sp_constraint[jrow].erase(it_other);
                }
            }
        }
    }
    // Remove emptry entries from the sp_constraint vector
    sp_constraint.erase(std::remove_if(sp_constraint.begin(),
                                       sp_constraint.end(),
                                       [](const std::map<unsigned int, double> &obj) { return obj.empty(); }),
                        sp_constraint.end());

    std::sort(sp_constraint.begin(), sp_constraint.end());
}


void rref_sparse2(const int ncols,
                  ConstraintSparseForm &sp_constraint,
                  const double tolerance)
{
    std::map<unsigned int, double> const_tmp2;
    int nrows = sp_constraint.size();
    int icol, irow, jrow;
    int jcol;
    double scaling_factor;
    double division_factor;

    int nrank = 0;

  /*  std::vector<std::set<unsigned int>> FirstColList;
    FirstColList.resize(ncols);*/

    // Stores pivot candidates of the input matrix
    //for (irow = 0; irow < nrows; ++irow) {
    //    FirstColList[sp_constraint[irow].begin()->first].insert(irow);
    //}

    unsigned int nsize_firstcol;

    int pivot;
    icol = 0;

    std::set<unsigned int>::iterator it_found;
    std::map<unsigned int, double>::iterator it_elem;

    for (irow = 0; irow < nrows; ++irow) {

        pivot = irow;

        while (true) {
            it_elem = sp_constraint[pivot].find(icol);
            if (it_elem != sp_constraint[pivot].end()) {
                if (std::abs(it_elem->second) > tolerance) {
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

        //for (jcol = icol; jcol < ncols; ++jcol) {
        //    it_found = std::find_if(FirstColList[jcol].begin(),
        //                            FirstColList[jcol].end(),
        //                            [irow](int x) { return x >= irow; });

        //    if (it_found != FirstColList[jcol].end()) {
        //        icol = jcol;
        //        break;
        //    }
        //}

        if (icol == ncols) break;

        if (std::abs(it_elem->second) > tolerance) ++nrank;
        //pivot = *it_found;

        if (pivot != irow) {

            /*   if (!sp_constraint[irow].empty()) {
                   FirstColList[sp_constraint[irow].begin()->first].erase(irow);
               }
               FirstColList[sp_constraint[pivot].begin()->first].erase(pivot);*/

            std::iter_swap(sp_constraint.begin() + irow,
                           sp_constraint.begin() + pivot);

            /*          if (!sp_constraint[pivot].empty()) {
                          FirstColList[sp_constraint[pivot].begin()->first].insert(pivot);
                      }
                      FirstColList[sp_constraint[irow].begin()->first].insert(irow);*/
        }

        division_factor = sp_constraint[irow].begin()->second;
        division_factor = 1.0 / division_factor;
        for (auto &it : sp_constraint[irow]) {
            it.second *= division_factor;
        }

        for (jrow = 0; jrow < nrows; ++jrow) {
            if (jrow == irow) continue;

            it_elem = sp_constraint[jrow].find(icol);
            if (it_elem == sp_constraint[jrow].end()) continue;
            scaling_factor = it_elem->second;

  /*          it_found = FirstColList[icol].find(jrow);
            if (it_found != FirstColList[icol].end()) {
                FirstColList[icol].erase(jrow);
            }*/

            // Subtract irow elements from jrow
            for (const auto &it_now : sp_constraint[irow]) {
                auto it_other = sp_constraint[jrow].find(it_now.first);
                if (it_other != sp_constraint[jrow].end()) {
                    it_other->second -= scaling_factor * it_now.second;
                } else {
                    sp_constraint[jrow][it_now.first] = -scaling_factor * it_now.second;
                }
            }
            for (auto it_other = sp_constraint[jrow].begin(); it_other != sp_constraint[jrow].end(); ++it_other) {
                if (std::abs(it_other->second) < tolerance) {
                    sp_constraint[jrow].erase(it_other);
                }
            }
            //if (!sp_constraint[jrow].empty()) {
            //    FirstColList[sp_constraint[jrow].begin()->first].insert(jrow);
            //}
        }
    }

    // Remove emptry entries from the sp_constraint vector
    sp_constraint.erase(std::remove_if(sp_constraint.begin(),
                                       sp_constraint.end(),
                                       [](const std::map<unsigned int, double> &obj) { return obj.empty(); }),
                        sp_constraint.end());

    // sp_constraint.erase(sp_constraint.begin() + nrank, sp_constraint.end());
    // sp_constraint.shrink_to_fit();
}
