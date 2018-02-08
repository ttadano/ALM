/*
 fitting.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include "fitting.h"
#include "files.h"
#include "error.h"
#include "memory.h"
#include "symmetry.h"
#include "system.h"
#include "fcs.h"
#include "interaction.h"
#include "timer.h"
#include "constants.h"
#include "constraint.h"
#include "mathfunctions.h"
#include <time.h>


using namespace ALM_NS;


Fitting::Fitting()
{
    set_default_variables();
}

Fitting::~Fitting()
{
    deallocate_variables();
}

void Fitting::set_default_variables()
{
    params = nullptr;
    u_in = nullptr;
    f_in = nullptr;
    ndata = 0;
    nstart = 1;
    nend = 0;
}

void Fitting::deallocate_variables()
{
    if (params) {
        deallocate(params);
    }
    if (u_in) {
        deallocate(u_in);
    }
    if (f_in) {
        deallocate(f_in);
    }
}

void Fitting::fitmain(ALM *alm)
{
    int i;
    int nat = alm->system->supercell.number_of_atoms;
    int natmin = alm->symmetry->nat_prim;
 //   int nstart = alm->fitting->nstart;
 //   int nend = alm->system->nend;
    int N, M, N_new;
    int maxorder = alm->interaction->maxorder;
    int P = alm->constraint->P;
    int ndata_used = nend - nstart + 1;

    double **u, **f;
    double **amat, *amat_1D, *fsum;
    double *fsum_orig;
    double *param_tmp;

    int nmulti = alm->symmetry->ntran;

    amat = nullptr;
    amat_1D = nullptr;
    fsum = nullptr;
    fsum_orig = nullptr;
    param_tmp = nullptr;

    alm->timer->start_clock("fitting");

    std::cout << " FITTING" << std::endl;
    std::cout << " =======" << std::endl << std::endl;

    std::cout << "  Reference files" << std::endl;
    std::cout << "   Displacement: " << alm->files->file_disp << std::endl;
    std::cout << "   Force       : " << alm->files->file_force << std::endl;
    std::cout << std::endl;

    std::cout << "  NSTART = " << nstart << "; NEND = " << nend << std::endl;
    std::cout << "  " << ndata_used << " entries will be used for fitting."
        << std::endl << std::endl;


    if (nmulti > 0) {
        allocate(u, ndata_used * nmulti, 3 * nat);
        allocate(f, ndata_used * nmulti, 3 * nat);
    } else {
        exit("fitmain", "nmulti has to be larger than 0.");
    }
    data_multiplier(u, f, nat, ndata_used, nmulti, alm->symmetry);

    N = 0;
    for (i = 0; i < maxorder; ++i) {
        N += alm->fcs->nequiv[i].size();
    }
    std::cout << "  Total Number of Parameters : "
        << N << std::endl << std::endl;

    // Calculate matrix elements for fitting

    M = 3 * natmin * ndata_used * nmulti;

    if (alm->constraint->constraint_algebraic) {
        N_new = 0;
        for (i = 0; i < maxorder; ++i) {
            N_new += alm->constraint->index_bimap[i].size();
        }
        std::cout << "  Total Number of Free Parameters : "
            << N_new << std::endl << std::endl;

        allocate(amat, M, N_new);
        allocate(fsum, M);
        allocate(fsum_orig, M);

        calc_matrix_elements_algebraic_constraint(M, N, N_new, nat, natmin, ndata_used,
                                                  nmulti, maxorder, u, f, amat, fsum,
                                                  fsum_orig, alm->symmetry, alm->fcs, alm->constraint);
    } else {
        allocate(amat, M, N);
        allocate(fsum, M);

        calc_matrix_elements(M, N, nat, natmin, ndata_used,
                             nmulti, maxorder, u, f, amat, fsum,
                             alm->symmetry, alm->fcs);
    }

    deallocate(u);
    deallocate(f);

    // Execute fitting

    allocate(param_tmp, N);


    // Fitting with singular value decomposition or QR-Decomposition

    if (alm->constraint->constraint_algebraic) {
        fit_algebraic_constraints(N_new, M, amat, fsum, param_tmp,
                                  fsum_orig, maxorder, alm->fcs, alm->constraint);

    } else if (alm->constraint->exist_constraint) {
        fit_with_constraints(N, M, P, amat, fsum, param_tmp,
                             alm->constraint->const_mat,
                             alm->constraint->const_rhs);
    } else {
        fit_without_constraints(N, M, amat, fsum, param_tmp);
    }


    // Copy force constants to public variable "params"
    if (params) {
        deallocate(params);
    }
    allocate(params, N);

    if (alm->constraint->constraint_algebraic) {

        for (i = 0; i < N; ++i) {
            params[i] = param_tmp[i];
        }
        deallocate(fsum_orig);

    } else {

        for (i = 0; i < N; ++i) params[i] = param_tmp[i];

    }

    if (amat) {
        deallocate(amat);
    }

    if (fsum) {
        deallocate(fsum);
    }

    if (param_tmp) {
        deallocate(param_tmp);
    }

    std::cout << std::endl;
    alm->timer->print_elapsed();
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    alm->timer->stop_clock("fitting");
}

void Fitting::set_displacement_and_force(const double * const *disp_in,
                                         const double * const *force_in,
                                         const int nat,
                                         const int ndata_used)
{
    if (u_in) {
        deallocate(u_in);
    }
    allocate(u_in, ndata_used, 3 * nat);

    if (f_in) {
        deallocate(f_in);
    }
    allocate(f_in, ndata_used, 3 * nat);

    for (int i = 0; i < ndata_used; i++) {
        for (int j = 0; j < 3 * nat; j++) {
            u_in[i][j] = disp_in[i][j];
            f_in[i][j] = force_in[i][j];
        }
    }
}

void Fitting::fit_without_constraints(int N,
                                      int M,
                                      double **amat,
                                      double *bvec,
                                      double *param_out)
{
    int i, j;
    unsigned long k;
    int nrhs = 1, nrank, INFO, LWORK;
    int LMIN, LMAX;
    double rcond = -1.0;
    double f_square = 0.0;
    double *WORK, *S, *amat_mod, *fsum2;

    std::cout << "  Entering fitting routine: SVD without constraints" << std::endl;

    LMIN = std::min<int>(M, N);
    LMAX = std::max<int>(M, N);

    LWORK = 3 * LMIN + std::max<int>(2 * LMIN, LMAX);
    LWORK = 2 * LWORK;

    allocate(WORK, LWORK);
    allocate(S, LMIN);

    // transpose matrix A
    allocate(amat_mod, M * N);
    allocate(fsum2, LMAX);

    k = 0;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            amat_mod[k++] = amat[i][j];
        }
    }
    for (i = 0; i < M; ++i) {
        fsum2[i] = bvec[i];
        f_square += std::pow(bvec[i], 2);
    }
    for (i = M; i < LMAX; ++i) fsum2[i] = 0.0;

    std::cout << "  SVD has started ... ";

    // Fitting with singular value decomposition
    dgelss_(&M, &N, &nrhs, amat_mod, &M, fsum2, &LMAX,
            S, &rcond, &nrank, WORK, &LWORK, &INFO);

    std::cout << "finished !" << std::endl << std::endl;

    std::cout << "  RANK of the matrix = " << nrank << std::endl;
    if (nrank < N)
        warn("fit_without_constraints",
             "Matrix is rank-deficient. Force constants could not be determined uniquely :(");

    if (nrank == N) {
        double f_residual = 0.0;
        for (i = N; i < M; ++i) {
            f_residual += std::pow(fsum2[i], 2);
        }
        std::cout << std::endl << "  Residual sum of squares for the solution: "
            << sqrt(f_residual) << std::endl;
        std::cout << "  Fitting error (%) : "
            << sqrt(f_residual / f_square) * 100.0 << std::endl;
    }

    for (i = 0; i < N; ++i) {
        param_out[i] = fsum2[i];
    }

    deallocate(WORK);
    deallocate(S);
    deallocate(fsum2);
    deallocate(amat_mod);
}

void Fitting::fit_with_constraints(int N,
                                   int M,
                                   int P,
                                   double **amat,
                                   double *bvec,
                                   double *param_out,
                                   double **cmat,
                                   double *dvec)
{
    int i, j;
    unsigned long k;
    int nrank;
    double f_square, f_residual;
    double *fsum2;
    double *mat_tmp;

    std::cout << "  Entering fitting routine: QRD with constraints" << std::endl;

    allocate(fsum2, M);
    allocate(mat_tmp, (M + P) * N);

    k = 0;

    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            mat_tmp[k++] = amat[i][j];
        }
        for (i = 0; i < P; ++i) {
            mat_tmp[k++] = cmat[i][j];
        }
    }

    nrank = rankQRD((M + P), N, mat_tmp, eps12);
    deallocate(mat_tmp);

    if (nrank != N) {
        std::cout << std::endl;
        std::cout << " **************************************************************************" << std::endl;
        std::cout << "  WARNING : rank deficient.                                                " << std::endl;
        std::cout << "  rank ( (A) ) ! = N            A: Fitting matrix     B: Constraint matrix " << std::endl;
        std::cout << "       ( (B) )                  N: The number of parameters                " << std::endl;
        std::cout << "  rank = " << nrank << " N = " << N << std::endl << std::endl;
        std::cout << "  This can cause a difficulty in solving the fitting problem properly      " << std::endl;
        std::cout << "  with DGGLSE, especially when the difference is large. Please check if    " << std::endl;
        std::cout << "  you obtain reliable force constants in the .fcs file.                    " << std::endl << std::endl;
        std::cout << "  You may need to reduce the cutoff radii and/or increase NDATA            " << std::endl;
        std::cout << "  by giving linearly-independent displacement patterns.                    " << std::endl;
        std::cout << " **************************************************************************" << std::endl;
        std::cout << std::endl;
    }

    f_square = 0.0;
    for (i = 0; i < M; ++i) {
        fsum2[i] = bvec[i];
        f_square += std::pow(bvec[i], 2);
    }
    std::cout << "  QR-Decomposition has started ...";

    double *amat_mod, *cmat_mod;
    allocate(amat_mod, M * N);
    allocate(cmat_mod, P * N);

    // transpose matrix A and C
    k = 0;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            amat_mod[k++] = amat[i][j];
        }
    }
    k = 0;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < P; ++i) {
            cmat_mod[k++] = cmat[i][j];
        }
    }

    // Fitting

    int LWORK = P + std::min<int>(M, N) + 10 * std::max<int>(M, N);
    int INFO;
    double *WORK, *x;
    allocate(WORK, LWORK);
    allocate(x, N);

    dgglse_(&M, &N, &P, amat_mod, &M, cmat_mod, &P,
            fsum2, dvec, x, WORK, &LWORK, &INFO);

    std::cout << " finished. " << std::endl;

    f_residual = 0.0;
    for (i = N - P; i < M; ++i) {
        f_residual += std::pow(fsum2[i], 2);
    }
    std::cout << std::endl << "  Residual sum of squares for the solution: "
        << sqrt(f_residual) << std::endl;
    std::cout << "  Fitting error (%) : "
        << std::sqrt(f_residual / f_square) * 100.0 << std::endl;

    // copy fcs to bvec

    for (i = 0; i < N; ++i) {
        param_out[i] = x[i];
    }

    deallocate(amat_mod);
    deallocate(cmat_mod);
    deallocate(WORK);
    deallocate(x);
    deallocate(fsum2);
}

void Fitting::fit_algebraic_constraints(int N,
                                        int M,
                                        double **amat,
                                        double *bvec,
                                        double *param_out,
                                        double *bvec_orig,
                                        const int maxorder,
                                        Fcs *fcs,
                                        Constraint *constraint)
{
    int i, j;
    unsigned long k;
    int nrhs = 1, nrank, INFO, LWORK;
    int LMIN, LMAX;
    double rcond = -1.0;
    double f_square = 0.0;
    double *WORK, *S, *amat_mod, *fsum2;

    std::cout << "  Entering fitting routine: SVD with constraints considered algebraically." << std::endl;

    LMIN = std::min<int>(M, N);
    LMAX = std::max<int>(M, N);

    LWORK = 3 * LMIN + std::max<int>(2 * LMIN, LMAX);
    LWORK = 2 * LWORK;

    allocate(WORK, LWORK);
    allocate(S, LMIN);

    // transpose matrix A
    allocate(amat_mod, M * N);
    allocate(fsum2, LMAX);

    k = 0;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            amat_mod[k++] = amat[i][j];
        }
    }
    for (i = 0; i < M; ++i) {
        fsum2[i] = bvec[i];
        f_square += std::pow(bvec_orig[i], 2);
    }
    for (i = M; i < LMAX; ++i) fsum2[i] = 0.0;

    std::cout << "  SVD has started ... ";

    // Fitting with singular value decomposition
    dgelss_(&M, &N, &nrhs, amat_mod, &M, fsum2, &LMAX,
            S, &rcond, &nrank, WORK, &LWORK, &INFO);

    std::cout << "finished !" << std::endl << std::endl;

    std::cout << "  RANK of the matrix = " << nrank << std::endl;
    if (nrank < N)
        warn("fit_without_constraints",
             "Matrix is rank-deficient. Force constants could not be determined uniquely :(");

    if (nrank == N) {
        double f_residual = 0.0;
        for (i = N; i < M; ++i) {
            f_residual += std::pow(fsum2[i], 2);
        }
        std::cout << std::endl;
        std::cout << "  Residual sum of squares for the solution: "
            << sqrt(f_residual) << std::endl;
        std::cout << "  Fitting error (%) : "
            << sqrt(f_residual / f_square) * 100.0 << std::endl;
    }

    int ishift = 0;
    int iparam = 0;
    double tmp;
    int inew, iold;

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < constraint->const_fix[i].size(); ++j) {
            param_out[constraint->const_fix[i][j].p_index_target + ishift]
                = constraint->const_fix[i][j].val_to_fix;
        }

        for (boost::bimap<int, int>::const_iterator it = constraint->index_bimap[i].begin();
             it != constraint->index_bimap[i].end(); ++it) {
            inew = (*it).left + iparam;
            iold = (*it).right + ishift;

            param_out[iold] = fsum2[inew];
        }

        for (j = 0; j < constraint->const_relate[i].size(); ++j) {
            tmp = 0.0;

            for (k = 0; k < constraint->const_relate[i][j].alpha.size(); ++k) {
                tmp += constraint->const_relate[i][j].alpha[k]
                    * param_out[constraint->const_relate[i][j].p_index_orig[k] + ishift];
            }
            param_out[constraint->const_relate[i][j].p_index_target + ishift] = -tmp;
        }

        ishift += fcs->nequiv[i].size();
        iparam += constraint->index_bimap[i].size();
    }

    deallocate(WORK);
    deallocate(S);
    deallocate(fsum2);
    deallocate(amat_mod);
}


void Fitting::calc_matrix_elements(const int M,
                                   const int N,
                                   const int nat,
                                   const int natmin,
                                   const int ndata_fit,
                                   const int nmulti,
                                   const int maxorder,
                                   double **u,
                                   double **f,
                                   double **amat,
                                   double *bvec,
                                   Symmetry *symmetry,
                                   Fcs *fcs)
{
    int i, j;
    int irow;
    int ncycle;

    std::cout << "  Calculation of matrix elements for direct fitting started ... ";
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            amat[i][j] = 0.0;
        }
        bvec[i] = 0.0;
    }

    ncycle = ndata_fit * nmulti;

#ifdef _OPENMP
#pragma omp parallel private(irow, i, j)
#endif
    {
        int *ind;
        int mm, order, iat, k;
        int im, idata, iparam;
        double amat_tmp;

        allocate(ind, maxorder + 1);

#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
        for (irow = 0; irow < ncycle; ++irow) {

            // generate r.h.s vector B
            for (i = 0; i < natmin; ++i) {
                iat = symmetry->map_p2s[i][0];
                for (j = 0; j < 3; ++j) {
                    im = 3 * i + j + 3 * natmin * irow;
                    bvec[im] = f[irow][3 * iat + j];
                }
            }

            // generate l.h.s. matrix A

            idata = 3 * natmin * irow;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                mm = 0;

                for (auto iter = fcs->nequiv[order].begin(); iter != fcs->nequiv[order].end(); ++iter) {
                    for (i = 0; i < *iter; ++i) {
                        ind[0] = fcs->fc_table[order][mm].elems[0];
                        k = idata + inprim_index(fcs->fc_table[order][mm].elems[0], symmetry);
                        amat_tmp = 1.0;
                        for (j = 1; j < order + 2; ++j) {
                            ind[j] = fcs->fc_table[order][mm].elems[j];
                            amat_tmp *= u[irow][fcs->fc_table[order][mm].elems[j]];
                        }
                        amat[k][iparam] -= gamma(order + 2, ind) * fcs->fc_table[order][mm].sign * amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }
        }

        deallocate(ind);

    }

    std::cout << "done!" << std::endl << std::endl;
}


void Fitting::calc_matrix_elements_algebraic_constraint(const int M,
                                                        const int N,
                                                        const int N_new,
                                                        const int nat,
                                                        const int natmin,
                                                        const int ndata_fit,
                                                        const int nmulti,
                                                        const int maxorder,
                                                        double **u,
                                                        double **f,
                                                        double **amat,
                                                        double *bvec,
                                                        double *bvec_orig,
                                                        Symmetry *symmetry,
                                                        Fcs *fcs,
                                                        Constraint *constraint)
{
    int i, j;
    int irow;
    int ncycle;

    std::cout << "  Calculation of matrix elements for direct fitting started ... ";

    ncycle = ndata_fit * nmulti;


#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N_new; ++j) {
            amat[i][j] = 0.0;
        }
        bvec[i] = 0.0;
        bvec_orig[i] = 0.0;
    }

#ifdef _OPENMP
#pragma omp parallel private(irow, i, j)
#endif
    {
        int *ind;
        int mm, order, iat, k;
        int im, idata, iparam;
        int ishift;
        int iold, inew;
        double amat_tmp;
        double **amat_orig;
        double **amat_mod;

        allocate(ind, maxorder + 1);
        allocate(amat_orig, 3 * natmin, N);
        allocate(amat_mod, 3 * natmin, N_new);

#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
        for (irow = 0; irow < ncycle; ++irow) {

            // generate r.h.s vector B
            for (i = 0; i < natmin; ++i) {
                iat = symmetry->map_p2s[i][0];
                for (j = 0; j < 3; ++j) {
                    im = 3 * i + j + 3 * natmin * irow;
                    bvec[im] = f[irow][3 * iat + j];
                    bvec_orig[im] = f[irow][3 * iat + j];
                }
            }

            for (i = 0; i < 3 * natmin; ++i) {
                for (j = 0; j < N; ++j) {
                    amat_orig[i][j] = 0.0;
                }
                for (j = 0; j < N_new; ++j) {
                    amat_mod[i][j] = 0.0;
                }
            }

            // generate l.h.s. matrix A

            idata = 3 * natmin * irow;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                mm = 0;

                for (auto iter = fcs->nequiv[order].begin(); iter != fcs->nequiv[order].end(); ++iter) {
                    for (i = 0; i < *iter; ++i) {
                        ind[0] = fcs->fc_table[order][mm].elems[0];
                        k = inprim_index(ind[0], symmetry);

                        amat_tmp = 1.0;
                        for (j = 1; j < order + 2; ++j) {
                            ind[j] = fcs->fc_table[order][mm].elems[j];
                            amat_tmp *= u[irow][fcs->fc_table[order][mm].elems[j]];
                        }
                        amat_orig[k][iparam] -= gamma(order + 2, ind) * fcs->fc_table[order][mm].sign * amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }

            ishift = 0;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                for (i = 0; i < constraint->const_fix[order].size(); ++i) {

                    for (j = 0; j < 3 * natmin; ++j) {
                        bvec[j + idata] -= constraint->const_fix[order][i].val_to_fix
                            * amat_orig[j][ishift + constraint->const_fix[order][i].p_index_target];
                    }
                }

                for (boost::bimap<int, int>::const_iterator it = constraint->index_bimap[order].begin();
                     it != constraint->index_bimap[order].end(); ++it) {
                    inew = (*it).left + iparam;
                    iold = (*it).right + ishift;

                    for (j = 0; j < 3 * natmin; ++j) {
                        amat_mod[j][inew] = amat_orig[j][iold];
                    }
                }

                for (i = 0; i < constraint->const_relate[order].size(); ++i) {

                    iold = constraint->const_relate[order][i].p_index_target + ishift;

                    for (j = 0; j < constraint->const_relate[order][i].alpha.size(); ++j) {

                        inew = constraint->index_bimap[order].right.at(
                                                                 constraint->const_relate[order][i].p_index_orig[j])
                            + iparam;
                        for (k = 0; k < 3 * natmin; ++k) {
                            amat_mod[k][inew] -= amat_orig[k][iold] * constraint->const_relate[order][i].alpha[j];
                        }
                    }
                }

                ishift += fcs->nequiv[order].size();
                iparam += constraint->index_bimap[order].size();
            }

            for (i = 0; i < 3 * natmin; ++i) {
                for (j = 0; j < N_new; ++j) {
                    amat[i + idata][j] = amat_mod[i][j];
                }
            }

        }

        deallocate(ind);
        deallocate(amat_orig);
        deallocate(amat_mod);
    }

    std::cout << "done!" << std::endl << std::endl;
}


void Fitting::data_multiplier(double **u,
                              double **f,
                              const int nat,
                              const int ndata_used,
                              const int nmulti,
                              Symmetry *symmetry)
{
    int i, j, k;
    int idata, itran, isym;
    int n_mapped;
    double u_rot[3], f_rot[3];

    // Multiply data
    idata = 0;
    for (i = 0; i < ndata_used; ++i) {
        for (itran = 0; itran < symmetry->ntran; ++itran) {
            for (j = 0; j < nat; ++j) {
                n_mapped = symmetry->map_sym[j][symmetry->symnum_tran[itran]];
                for (k = 0; k < 3; ++k) {
                    u[idata][3 * n_mapped + k] = u_in[i][3 * j + k];
                    f[idata][3 * n_mapped + k] = f_in[i][3 * j + k];
                }
            }
            ++idata;
        }
    }
}


int Fitting::inprim_index(const int n, Symmetry *symmetry)
{
    int in;
    int atmn = n / 3;
    int crdn = n % 3;

    for (int i = 0; i < symmetry->nat_prim; ++i) {
        if (symmetry->map_p2s[i][0] == atmn) {
            in = 3 * i + crdn;
            break;
        }
    }
    return in;
}

double Fitting::gamma(const int n, const int *arr)
{
    int *arr_tmp, *nsame;
    int i;
    int ind_front, nsame_to_front;

    allocate(arr_tmp, n);
    allocate(nsame, n);

    for (i = 0; i < n; ++i) {
        arr_tmp[i] = arr[i];
        nsame[i] = 0;
    }

    ind_front = arr[0];
    nsame_to_front = 1;

    insort(n, arr_tmp);

    int nuniq = 1;
    int iuniq = 0;

    nsame[0] = 1;

    for (i = 1; i < n; ++i) {
        if (arr_tmp[i] == arr_tmp[i - 1]) {
            ++nsame[iuniq];
        } else {
            ++nsame[++iuniq];
            ++nuniq;
        }

        if (arr[i] == ind_front) ++nsame_to_front;
    }

    int denom = 1;

    for (i = 0; i < nuniq; ++i) {
        denom *= factorial(nsame[i]);
    }

    deallocate(arr_tmp);
    deallocate(nsame);

    return static_cast<double>(nsame_to_front) / static_cast<double>(denom);
}

int Fitting::factorial(const int n)
{
    if (n == 1 || n == 0) {
        return 1;
    } else {
        return n * factorial(n - 1);
    }
}


int Fitting::rankQRD(const int m,
                     const int n,
                     double *mat,
                     const double tolerance)
{
    // Return the rank of matrix mat revealed by the column pivoting QR decomposition
    // The matrix mat is destroyed.

    int m_ = m;
    int n_ = n;

    int LDA = m_;

    int LWORK = 10 * n_;
    int INFO;
    int *JPVT;
    double *WORK, *TAU;

    int nmin = std::min<int>(m_, n_);

    allocate(JPVT, n_);
    allocate(WORK, LWORK);
    allocate(TAU, nmin);

    for (int i = 0; i < n_; ++i) JPVT[i] = 0;

    dgeqp3_(&m_, &n_, mat, &LDA, JPVT, TAU, WORK, &LWORK, &INFO);

    deallocate(JPVT);
    deallocate(WORK);
    deallocate(TAU);

    if (std::abs(mat[0]) < eps) return 0;

    double **mat_tmp;
    allocate(mat_tmp, m_, n_);

    unsigned long k = 0;

    for (int j = 0; j < n_; ++j) {
        for (int i = 0; i < m_; ++i) {
            mat_tmp[i][j] = mat[k++];
        }
    }

    int nrank = 0;
    for (int i = 0; i < nmin; ++i) {
        if (std::abs(mat_tmp[i][i]) > tolerance * std::abs(mat[0])) ++nrank;
    }

    deallocate(mat_tmp);

    return nrank;
}

int Fitting::rankSVD(const int m,
                     const int n,
                     double *mat,
                     const double tolerance)
{
    int i;
    int m_ = m;
    int n_ = n;

    int LWORK = 10 * m;
    int INFO;
    int *IWORK;
    int ldu = 1, ldvt = 1;
    double *s, *WORK;
    double u[1], vt[1];

    int nmin = std::min<int>(m, n);

    allocate(IWORK, 8 * nmin);
    allocate(WORK, LWORK);
    allocate(s, nmin);

    char mode[] = "N";

    dgesdd_(mode, &m_, &n_, mat, &m_, s, u, &ldu, vt, &ldvt,
            WORK, &LWORK, IWORK, &INFO);

    int rank = 0;
    for (i = 0; i < nmin; ++i) {
        if (s[i] > s[0] * tolerance) ++rank;
    }

    deallocate(WORK);
    deallocate(IWORK);
    deallocate(s);

    return rank;
}

int Fitting::rankSVD2(const int m_in,
                      const int n_in,
                      double **mat,
                      const double tolerance)
{
    // Reveal the rank of matrix mat without destroying the matrix elements

    int i, j, k;
    double *arr;

    int m = m_in;
    int n = n_in;

    allocate(arr, m * n);

    k = 0;

    for (j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i) {
            arr[k++] = mat[i][j];
        }
    }

    int LWORK = 10 * m;
    int INFO;
    int *IWORK;
    int ldu = 1, ldvt = 1;
    double *s, *WORK;
    double u[1], vt[1];

    int nmin = std::min<int>(m, n);

    allocate(IWORK, 8 * nmin);
    allocate(WORK, LWORK);
    allocate(s, nmin);

    char mode[] = "N";

    dgesdd_(mode, &m, &n, arr, &m, s, u, &ldu, vt, &ldvt,
            WORK, &LWORK, IWORK, &INFO);

    int rank = 0;
    for (i = 0; i < nmin; ++i) {
        if (s[i] > s[0] * tolerance) ++rank;
    }

    deallocate(IWORK);
    deallocate(WORK);
    deallocate(s);
    deallocate(arr);

    return rank;
}
