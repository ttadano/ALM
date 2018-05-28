/*
lasso.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "lasso.h"
#include "fitting.h"
#include "memory.h"
#include "system.h"
#include "symmetry.h"
#include "files.h"
#include "interaction.h"
#include "input_parser.h"
#include "fcs.h"
#include "constraint.h"
#include "timer.h"
#include "error.h"
#include <iostream>
#include <iomanip>
#include <random>
#include "constants.h"
#include <Eigen/Dense>


using namespace ALM_NS;

Lasso::Lasso()
{
    set_default_variables();
}

Lasso::~Lasso()
{
}

void Lasso::set_default_variables()
{
    disp_norm = 1.0;
    l1_alpha = 1.0;
    l2_lambda = 10.0;
    lasso_tol = 1.0e-7;
    maxiter = 100000;
    maxiter_cg = 5;
    lasso_cv = 0;
    lasso_cvset = 10;
    output_frequency = 1000;
    lasso_zero_thr = 1.0e-50;
    l1_alpha_min = 1.0e-3;
    l1_alpha_max = 1.0;
    num_l1_alpha = 100;
    lasso_pcg = 0;
    lasso_algo = 0;
    standardize = 1;
}

void Lasso::lasso_main(ALM *alm)
{
    int i, j, k;
    const int nat = alm->system->supercell.number_of_atoms;
    const auto natmin = alm->symmetry->nat_prim;
    const int maxorder = alm->interaction->maxorder;
    const int ndata = alm->fitting->ndata;
    const int nstart = alm->fitting->nstart;
    const int nend = alm->fitting->nend;
    const int skip_s = alm->fitting->skip_s;
    const int skip_e = alm->fitting->skip_e;
    const int ntran = alm->symmetry->ntran;
    const int ndata_used = nend - nstart + 1 - skip_e + skip_s;

    double scale_factor;
    double **u, **f;
    double **amat, *fsum;
    double *bvec_breg, *dvec_breg;
    // for LASSO validation
    double **u_test, **f_test;
    double **amat_test, *fsum_test;
    double fnorm, fnorm_test;
    int ndata_used_test = alm->lasso->nend_test - alm->lasso->nstart_test + 1;


    int N = 0;
    int N_new = 0;
    for (auto i = 0; i < maxorder; ++i) {
        N += alm->fcs->nequiv[i].size();
        N_new += alm->constraint->index_bimap[i].size();
    }

    int M = 3 * natmin * static_cast<long>(ndata_used) * ntran;
    int M_test = 3 * natmin * ndata_used_test * ntran;

    if (alm->verbosity > 0) {
        std::cout << " LASSO" << std::endl;
        std::cout << " =====" << std::endl << std::endl;

        std::cout << "  Reference files" << std::endl;
        std::cout << "   Displacement: " << alm->files->file_disp << std::endl;
        std::cout << "   Force       : " << alm->files->file_force << std::endl;
        std::cout << std::endl;

        std::cout << "  NSTART = " << nstart << "; NEND = " << nend;
        if (skip_s < skip_e) std::cout << ": SKIP = " << skip_s + 1 << "-" << skip_e;
        std::cout << std::endl;
        std::cout << "  " << ndata_used
            << " entries will be used for lasso." << std::endl << std::endl;

        std::cout << "  Validation test files" << std::endl;
        std::cout << "   Displacement: " << alm->lasso->dfile_test << std::endl;
        std::cout << "   Force       : " << alm->lasso->ffile_test << std::endl;
        std::cout << std::endl;

        std::cout << "  NSTART = " << nstart_test << "; NEND = " << nend_test << std::endl;
        std::cout << "  " << ndata_used_test
            << " entries will be used for lasso validation." << std::endl << std::endl;

        std::cout << "  Total Number of Parameters : " << N << std::endl;
        std::cout << "  Total Number of Free Parameters : " << N_new << std::endl << std::endl;
    }

    // Parse displacement and force
    allocate(u, ndata_used, 3 * nat);
    allocate(f, ndata_used, 3 * nat);
    allocate(u_test, ndata_used_test, 3 * nat);
    allocate(f_test, ndata_used_test, 3 * nat);

    InputParser *input_parser = new InputParser();

    input_parser->parse_displacement_and_force_files(u,
                                                     f,
                                                     nat,
                                                     ndata,
                                                     nstart,
                                                     nend,
                                                     skip_s,
                                                     skip_e,
                                                     alm->files->file_disp,
                                                     alm->files->file_force);

    input_parser->parse_displacement_and_force_files(u_test,
                                                     f_test,
                                                     nat,
                                                     ndata_used_test,
                                                     nstart_test,
                                                     nend_test,
                                                     0,
                                                     0,
                                                     dfile_test,
                                                     ffile_test);

    delete input_parser;

    // Scale displacements

    const double inv_dnorm = 1.0 / disp_norm;
    for (i = 0; i < ndata_used; ++i) {
        for (j = 0; j < 3 * nat; ++j) {
            u[i][j] *= inv_dnorm;
        }
    }

    for (i = 0; i < ndata_used_test; ++i) {
        for (j = 0; j < 3 * nat; ++j) {
            u_test[i][j] *= inv_dnorm;
        }
    }

    // Scale force constants
    for (i = 0; i < maxorder; ++i) {
        scale_factor = std::pow(disp_norm, i + 1);
        for (j = 0; j < alm->constraint->const_fix[i].size(); ++j) {
            alm->constraint->const_fix[i][j].val_to_fix *= scale_factor;
        }
    }


    unsigned long nrows = 3 * static_cast<long>(natmin)
        * static_cast<long>(ndata_used)
        * static_cast<long>(ntran);

    const unsigned long ncols = static_cast<long>(N_new);

    std::vector<double> amat_1D;
    std::vector<double> bvec;

    amat_1D.resize(nrows * ncols, 0.0);
    bvec.resize(nrows, 0.0);

    alm->fitting->set_displacement_and_force(u, f, nat, ndata_used);

    alm->fitting->get_matrix_elements_algebraic_constraint(maxorder,
                                                           ndata_used,
                                                           &amat_1D[0],
                                                           &bvec[0],
                                                           fnorm,
                                                           alm->symmetry,
                                                           alm->fcs,
                                                           alm->constraint);

    deallocate(u);
    deallocate(f);
    allocate(amat, M, N_new);
    allocate(fsum, M);

    unsigned long m = 0;
    for (auto icol = 0; icol < N_new; ++icol) {
        for (auto irow = 0; irow < M; ++irow) {
            amat[irow][icol] = amat_1D[m++];
        }
    }
    for (auto irow = 0; irow < M; ++irow) {
        fsum[irow] = bvec[irow];
    }

    nrows = 3 * static_cast<long>(natmin)
        * static_cast<long>(ndata_used_test)
        * static_cast<long>(ntran);

    amat_1D.resize(nrows * ncols, 0.0);
    bvec.resize(nrows, 0.0);

    alm->fitting->set_displacement_and_force(u_test, f_test, nat, ndata_used_test);

    alm->fitting->get_matrix_elements_algebraic_constraint(maxorder,
                                                           ndata_used_test,
                                                           &amat_1D[0],
                                                           &bvec[0],
                                                           fnorm_test,
                                                           alm->symmetry,
                                                           alm->fcs,
                                                           alm->constraint);

    deallocate(u_test);
    deallocate(f_test);
    allocate(amat_test, M_test, N_new);
    allocate(fsum_test, M_test);

    m = 0;
    for (auto icol = 0; icol < N_new; ++icol) {
        for (auto irow = 0; irow < M_test; ++irow) {
            amat_test[irow][icol] = amat_1D[m++];
        }
    }
    for (auto irow = 0; irow < M_test; ++irow) {
        fsum_test[irow] = bvec[irow];
    }

    amat_1D.clear();
    bvec.clear();


    // Scale back force constants

    for (i = 0; i < maxorder; ++i) {
        scale_factor = 1.0 / std::pow(disp_norm, i + 1);
        for (j = 0; j < alm->constraint->const_fix[i].size(); ++j) {
            alm->constraint->const_fix[i][j].val_to_fix *= scale_factor;
        }
    }

    // Start Lasso optimization

    std::vector<double> param(N_new);

    double *factor_std;
    bool *has_prod;

    Eigen::MatrixXd A, Prod;
    Eigen::VectorXd b, C, grad, x;
    Eigen::VectorXd scale_beta;

    if (lasso_algo == 0) {

        // Coordinate descent

        A.resize(M, N_new);
        Prod.resize(N_new, N_new);
        b.resize(M);
        C.resize(N_new);
        grad.resize(N_new);
        x.resize(N_new);
        scale_beta.resize(N_new);

        allocate(has_prod, N_new);
        allocate(factor_std, N_new);

        for (i = 0; i < N_new; ++i) {
            param[i] = 0.0;
            x[i] = 0.0;
            has_prod[i] = false;
        }

        for (i = 0; i < M; ++i) {
            for (j = 0; j < N_new; ++j) {
                A(i, j) = amat[i][j];
            }
            b(i) = fsum[i];
        }

        // Standardize if necessary

        double Minv = 1.0 / static_cast<double>(M);

        if (standardize) {
            double sum1, sum2, tmp;

            std::cout << " STANDARDIZE = 1 : Standardization will be performed for matrix A and vector b." << std::endl;
            std::cout << "                   The LASSO_DNORM-tag will be neglected." << std::endl;
            for (j = 0; j < N_new; ++j) {
                sum1 = A.col(j).sum() * Minv;
                sum2 = A.col(j).dot(A.col(j)) * Minv;

                for (i = 0; i < M; ++i) {
                    A(i, j) = (A(i, j) - sum1) / std::sqrt(sum2 - sum1 * sum1);
                }
                factor_std[j] = 1.0 / std::sqrt(sum2 - sum1 * sum1);
                scale_beta(j) = 1.0;
            }
            // Is it necesary?
            tmp = b.sum() * Minv;
            for (i = 0; i < M; ++i) b(i) -= tmp;
        } else {
            double sum2;
            std::cout << " STANDARDIZE = 0 : No standardization of matrix A and vector b." << std::endl;
            std::cout << "                   Columns of matrix A will be scaled by the LASSO_DNORM value." << std::endl;
            for (j = 0; j < N_new; ++j) {
                factor_std[j] = 1.0;
                sum2 = A.col(j).dot(A.col(j)) * Minv;
                scale_beta(j) = 1.0 / sum2;
            }
        }

        C = A.transpose() * b;
        auto lambda_max = 0.0;
        for (i = 0; i < N_new; ++i) {
            lambda_max = std::max<double>(lambda_max, std::abs(C(i)));
        }
        lambda_max /= static_cast<double>(M);
        std::cout << std::endl;
        std::cout << " Recommended LASSO_MAXALPHA = " << lambda_max << std::endl << std::endl;
        grad = C;
        for (i = 0; i < N_new; ++i) {
            for (j = 0; j < N_new; ++j) {
                Prod(i, j) = 0.0;
            }
        }

    } else if (lasso_algo == 1) {

        // Split bregman

        allocate(bvec_breg, N_new);
        allocate(dvec_breg, N_new);

        for (i = 0; i < N_new; ++i) {
            param[i] = 0.0;
            bvec_breg[i] = 0.0;
            dvec_breg[i] = 0.0;
        }
    }

    if (lasso_cv == 1) {

        // Cross-validation mode

        std::cout << "  Lasso validation with the following parameters:" << std::endl;
        std::cout << "   LASSO_MINALPHA = " << std::setw(15) << l1_alpha_min;
        std::cout << " LASSO_MAXALPHA = " << std::setw(15) << l1_alpha_max << std::endl;
        std::cout << "   LASSO_NALPHA = " << std::setw(5) << num_l1_alpha << std::endl;
        std::cout << "   LASSO_TOL = " << std::setw(15) << lasso_tol << std::endl;
        std::cout << "   LASSO_MAXITER = " << std::setw(5) << maxiter << std::endl;
        std::cout << "   LASSO_DBASIS = " << std::setw(15) << disp_norm << std::endl;
        if (lasso_algo == 1) {
            std::cout << "   LASSO_LAMBDA (L2) = " << std::setw(15) << l2_lambda << std::endl;
            std::cout << "   LASSO_ZERO_THR = " << std::setw(15) << lasso_zero_thr << std::endl;
        }
        std::cout << std::endl;

        std::ofstream ofs_cv;

        std::string file_cv = alm->files->job_title + ".lasso_cv";
        ofs_cv.open(file_cv.c_str(), std::ios::out);

        if (lasso_algo == 0) {
            ofs_cv << "# Algorithm : Coordinate descent" << std::endl;
            ofs_cv << "# LASSO_DBASIS = " << std::setw(15) << disp_norm << std::endl;
            ofs_cv << "# LASSO_TOL = " << std::setw(15) << lasso_tol << std::endl;
            ofs_cv << "# L1 ALPHA, Fitting error, Validation error, Num. zero IFCs (2nd, 3rd, ...) " << std::endl;
        } else if (lasso_algo == 1) {
            ofs_cv << "# Algorithm : Split-Bregman iteration" << std::endl;
            ofs_cv << "# LASSO_LAMBDA (L2) = " << std::setw(15) << l2_lambda << std::endl;
            ofs_cv << "# LASSO_DBASIS = " << std::setw(15) << disp_norm << std::endl;
            ofs_cv << "# LASSO_TOL = " << std::setw(15) << lasso_tol << std::endl;
            ofs_cv << "# L1 ALPHA, Fitting error, Validation error, Num. zero IFCs (2nd, 3rd, ...) " << std::endl;
        }

        int initialize_mode;
        double res1, res2;

        std::vector<int> nzero_lasso(maxorder);

        for (int ialpha = 0; ialpha <= num_l1_alpha; ++ialpha) {

            l1_alpha = l1_alpha_min * std::pow(l1_alpha_max / l1_alpha_min,
                                               static_cast<double>(num_l1_alpha - ialpha) / static_cast<double>(
                                                   num_l1_alpha));

            std::cout << "-----------------------------------------------------------------" << std::endl;
            std::cout << "  L1_ALPHA = " << std::setw(15) << l1_alpha << std::endl;

            ofs_cv << std::setw(15) << l1_alpha;

            if (ialpha == 0) {
                initialize_mode = 0;
            } else {
                initialize_mode = 1;
            }

            if (lasso_algo == 0) {
                // Coordinate Descent method
                coordinate_descent(M, N_new, l1_alpha, lasso_tol, initialize_mode, maxiter,
                                   x, A, b, C, has_prod, Prod, grad, fnorm, output_frequency,
                                   scale_beta, standardize);

                for (i = 0; i < N_new; ++i) param[i] = x[i] * factor_std[i];


            } else if (lasso_algo == 1) {
                // Split-Bregman
                split_bregman_minimization(M, N_new, l1_alpha, l2_lambda, lasso_tol, maxiter,
                                           amat, fsum, fnorm, &param[0], bvec_breg, dvec_breg,
                                           initialize_mode, output_frequency);
            }

            calculate_residual(M, N_new, amat, &param[0], fsum, fnorm, res1);
            calculate_residual(M_test, N_new, amat_test, &param[0], fsum_test, fnorm_test, res2);

            // Count the number of zero parameters
            int iparam = 0;

            for (i = 0; i < maxorder; ++i) {
                nzero_lasso[i] = 0;
                for (const auto &it : alm->constraint->index_bimap[i]) {
                    int inew = it.left + iparam;
                    if (std::abs(param[inew]) < eps) ++nzero_lasso[i];

                }
                iparam += alm->constraint->index_bimap[i].size();
            }

            ofs_cv << std::setw(15) << std::sqrt(res1);
            ofs_cv << std::setw(15) << std::sqrt(res2);
            for (i = 0; i < maxorder; ++i) {
                ofs_cv << std::setw(10) << nzero_lasso[i];
            }
            ofs_cv << std::endl;
        }

        ofs_cv.close();

    } else if (lasso_cv == 0) {

        double res1;
        std::vector<int> nzero_lasso(maxorder);
        bool find_sparse_representation = true;

        std::cout << "  Lasso minimization with the following parameters:" << std::endl;
        std::cout << "   LASSO_ALPHA  (L1) = " << std::setw(15) << l1_alpha << std::endl;
        std::cout << "   LASSO_TOL = " << std::setw(15) << lasso_tol << std::endl;
        std::cout << "   LASSO_MAXITER = " << std::setw(5) << maxiter << std::endl;
        std::cout << "   LASSO_DBASIS = " << std::setw(15) << disp_norm << std::endl;
        if (lasso_algo == 1) {
            std::cout << "   LASSO_LAMBDA (L2) = " << std::setw(15) << l2_lambda << std::endl;
            std::cout << "   LASSO_ZERO_THR = " << std::setw(15) << lasso_zero_thr << std::endl;
        }
        std::cout << std::endl;

        if (lasso_algo == 0) {
            // Coordinate Descent Method
            find_sparse_representation = false;

            coordinate_descent(M, N_new, l1_alpha, lasso_tol, 0, maxiter,
                               x, A, b, C, has_prod, Prod, grad, fnorm * fnorm,
                               output_frequency, scale_beta, standardize);

            for (i = 0; i < N_new; ++i) param[i] = x[i] * factor_std[i];

        } else if (lasso_algo == 1) {
            // Split-Bregman Method
            split_bregman_minimization(M, N_new, l1_alpha, l2_lambda, lasso_tol, maxiter,
                                       amat, fsum, fnorm * fnorm, &param[0],
                                       bvec_breg, dvec_breg, 0, output_frequency);
        }


        calculate_residual(M, N_new, amat, &param[0], fsum, fnorm, res1);

        // Count the number of zero parameters
        int iparam, inew;
        iparam = 0;

        for (i = 0; i < maxorder; ++i) {
            nzero_lasso[i] = 0;
            for (boost::bimap<int, int>::const_iterator it = alm->constraint->index_bimap[i].begin();
                 it != alm->constraint->index_bimap[i].end(); ++it) {
                inew = (*it).left + iparam;
                if (std::abs(param[inew]) < eps) ++nzero_lasso[i];

            }
            iparam += alm->constraint->index_bimap[i].size();
        }

        std::cout << "  RESIDUAL (%): " << std::sqrt(res1) * 100.0 << std::endl;
        for (int order = 0; order < maxorder; ++order) {
            std::cout << "  Number of non-zero " << std::setw(9) << alm->interaction->str_order[order] << " FCs : "
                << alm->constraint->index_bimap[order].size() - nzero_lasso[order] << std::endl;
        }
        std::cout << std::endl;

        if (find_sparse_representation) {

            // Calculate prefactor for atomic forces
            std::vector<double> prefactor_force(N_new);

            get_prefactor_force(maxorder,
                                alm->fcs,
                                alm->constraint,
                                alm->fitting,
                                prefactor_force);

            std::vector<double> param_copy(N_new);
            double tolerance_tmp;

            for (i = 0; i < N_new; ++i) {
                param_copy[i] = param[i];
            }

            scale_factor = std::pow(1.0e+21, 0.01);

            std::cout << " LASSO_ZERO_THR, Number of zero parameters, fitting error (%) " << std::endl;

            for (int itol = 0; itol <= 100; ++itol) {
                tolerance_tmp = 1.0e-20 * std::pow(scale_factor, static_cast<double>(itol));

                std::cout << std::setw(15) << tolerance_tmp;

                iparam = 0;

                for (i = 0; i < maxorder; ++i) {
                    nzero_lasso[i] = 0;

                    for (const auto &it : alm->constraint->index_bimap[i]) {
                        inew = it.left + iparam;

                        if (std::abs(param_copy[inew]) * prefactor_force[inew] < tolerance_tmp) {
                            ++nzero_lasso[i];
                            param_copy[inew] = 0.0;
                        }
                    }
                    iparam += alm->constraint->index_bimap[i].size();

                    std::cout << std::setw(10) << nzero_lasso[i];
                }

                calculate_residual(M, N_new, amat, &param_copy[0], fsum, fnorm * fnorm, res1);
                std::cout << std::setw(15) << std::sqrt(res1) * 100.0 << std::endl;
            }
        }
    }

    if (lasso_algo == 0) {
        deallocate(has_prod);
        deallocate(factor_std);
    } else if (lasso_algo == 1) {
        deallocate(bvec_breg);
        deallocate(dvec_breg);
    }

    k = 0;
    for (i = 0; i < maxorder; ++i) {
        scale_factor = 1.0 / std::pow(disp_norm, i + 1);

        for (j = 0; j < alm->constraint->index_bimap[i].size(); ++j) {
            param[k] *= scale_factor;
            ++k;
        }
    }

    alm->fitting->set_fcs_values(maxorder,
                                 &param[0],
                                 alm->fcs->nequiv,
                                 alm->constraint);

    deallocate(amat);
    deallocate(fsum);
    deallocate(amat_test);
    deallocate(fsum_test);

    alm->timer->print_elapsed();
    std::cout << " --------------------------------------------------------------" << std::endl;
}


void Lasso::split_bregman_minimization(const int M,
                                       const int N,
                                       const double alpha,
                                       const double lambda,
                                       const double tolerance,
                                       const int maxiter,
                                       double **Amat,
                                       double *fvec,
                                       const double f2norm,
                                       double *param,
                                       double *bvec,
                                       double *dvec,
                                       const int cv_mode,
                                       const int nfreq)
{
    // Optimization of L1-penalized problem using split-Bregman algorithm

    int iter_lasso;
    int i, j;

    double tmp, diff;
    double invlambda = 1.0 / lambda;
    double al = alpha * lambda;
    bool do_print_log;
    int preconditioner = lasso_pcg;
    int miniter = 5;

    std::cout << "  Start LASSO minimization with the split-Bregman algorithm" << std::endl;

#ifdef _USE_EIGEN_DISABLED
    using namespace Eigen;

    MatrixXd Qmat(N, N);
    MatrixXd Amat2(M, N);
    MatrixXd mat_precondition(N, N);
    VectorXd fvec2(M), bvec_CG2(N), bvec_CG_update2(N);
    VectorXd bvec2(N), dvec2(N), param2(N);
    VectorXd res(N);
    VectorXd vec_d(N);

    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            Amat2(i, j) = Amat[i][j];
        }
        fvec2(i) = fvec[i];
    }

    // Prepare Qmat = A^T A, which is positive semidefinite
    // When alpha is nonzero, Qmat becomes positive definite.
    Qmat = Amat2.transpose() * Amat2;
    for (i = 0; i < N; ++i) {
        Qmat(i, i) += alpha * alpha * lambda;
    }

    bvec_CG2 = Amat2.transpose() * fvec2;

    if (preconditioner == 1) {
    //
    // Incomplete Cholesky factorization without fill-in, known as IC(0).
    //
        std::cout << "    LASSO_PCG = 1: Prepare preconditioning matrix ...";
        incomplete_cholesky_factorization(N, Qmat, mat_precondition, vec_d);
        std::cout << " done." << std::endl << std::endl;
    }

    if (cv_mode == 0) {
        std::cout << "    Start with b = 0 and d = 0" << std::endl;

        for (i = 0; i < N; ++i) {
            bvec2(i) = 0.0;
            dvec2(i) = 0.0;
            param2(i) = 0.0;
        }

    } else if (cv_mode == 1) {

        std::cout << "    Start from b, d, and x of the previous run" << std::endl;
        for (i = 0; i < N; ++i) {
            bvec2(i) = bvec[i];
            dvec2(i) = dvec[i];
            param2(i) = param[i];
        }
    } else if (cv_mode == 2) {

        std::cout << "    Start with b = 0 and d = 0" << std::endl;
        std::cout << "    Parameter x is initialized by solving the LS equation" << std::endl;

        for (i = 0; i < N; ++i) {
            bvec2(i) = 0.0;
            dvec2(i) = 0.0;
            param2(i) = 0.0;
        }

        bvec_CG_update2 = bvec_CG2 + alpha * lambda * (dvec2 - bvec2);
        minimize_quadratic_CG(N, Qmat, bvec_CG_update2, param2, 20 * N, true,
                              mat_precondition, vec_d, preconditioner);
    }

    for (iter_lasso = 0; iter_lasso < maxiter; ++iter_lasso) {

        do_print_log = !((iter_lasso + 1) % nfreq);

        if (do_print_log) {
            std::cout << "   SPLIT BREGMAN : " << std::setw(5) << iter_lasso + 1 << std::endl;
            std::cout << "    Update parameter vector using conjugate gradient" << std::endl;
        }

    // Update bvec_CG for CG
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < N; ++i) {
            bvec_CG_update2(i) = bvec_CG2(i) + al * (dvec2(i) - bvec2(i));
        }

    // Update parameters via CG minimization

        minimize_quadratic_CG(N, Qmat, bvec_CG_update2, param2, maxiter_cg, do_print_log,
                              mat_precondition, vec_d, preconditioner);

    // Update dvec using the shrink operator
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < N; ++i) {
            dvec2(i) = shrink(alpha * param2(i) + bvec2(i), invlambda);
        }

    // Update bvec
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < N; ++i) {
            bvec2(i) += alpha * param2(i) - dvec2(i);
        }

        diff = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:diff)
#endif
        for (i = 0; i < N; ++i) {
            diff += std::pow(param[i] - param2(i), 2);
        }

        if (do_print_log) {

            param2norm = param2.dot(param2);

            std::cout << std::endl;
            std::cout << "    1: ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << std::sqrt(diff / static_cast<double>(N))
                << std::setw(15) << std::sqrt(diff / param2norm) << std::endl;

            tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
            for (i = 0; i < N; ++i) {
                tmp += std::abs(param2(i));
            }
            std::cout << "    2: ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
            res = Amat2 * param2 - fvec2;
            tmp = res.dot(res);
            std::cout << "    3: ||Au_{k}-f||_2          = " << std::setw(15) << std::sqrt(tmp) << std::setw(15) << std::sqrt(tmp / f2norm) << std::endl;
            res = dvec2 - alpha * param2;
            std::cout << "    4: ||d_{k}-alpha*u_{k}||_2 = " << std::setw(15) << std::sqrt(res.dot(res)) << std::endl;

            tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
            for (i = 0; i < N; ++i) {
                tmp += std::abs(dvec2(i));
            }
            std::cout << "    5: ||d_{k}||_1             = " << std::setw(15) << tmp << std::endl;
            std::cout << "    6: ||b_{k}||_2             = " << std::setw(15) << std::sqrt(bvec2.dot(bvec2)) << std::endl;
            std::cout << std::endl;
        }

        for (i = 0; i < N; ++i) {
            param[i] = param2(i);
            bvec[i] = bvec2(i);
            dvec[i] = dvec2(i);
        }


        if (std::sqrt(diff / static_cast<double>(N)) < tolerance && iter_lasso >= miniter) {

            std::cout << "   SPLIT BREGMAN : " << std::setw(5) << iter_lasso + 1 << std::endl;
            std::cout << "    Update parameter vector using conjugate gradient" << std::endl;

            std::cout << std::endl;
            std::cout << "    1': ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << std::sqrt(diff / static_cast<double>(N))
                << std::setw(15) << std::sqrt(diff / param2.dot(param2)) << std::endl;

            tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
            for (i = 0; i < N; ++i) {
                tmp += std::abs(param2(i));
            }
            std::cout << "    2': ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
            res = Amat2 * param2 - fvec2;
            tmp = res.dot(res);
            std::cout << "    3': ||Au_{k}-f||_2          = " << std::setw(15) << std::sqrt(tmp) << std::setw(15) << std::sqrt(tmp / f2norm) << std::endl;
            res = dvec2 - alpha * param2;
            std::cout << "    4': ||d_{k}-alpha*u_{k}||_2 = " << std::setw(15) << std::sqrt(res.dot(res)) << std::endl;

            tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
            for (i = 0; i < N; ++i) {
                tmp += std::abs(dvec2(i));
            }
            std::cout << "    5': ||d_{k}||_1             = " << std::setw(15) << tmp << std::endl;
            std::cout << "    6': ||b_{k}||_2             = " << std::setw(15) << std::sqrt(bvec2.dot(bvec2)) << std::endl;
            std::cout << std::endl;

            for (i = 0; i < N; ++i) {
                param[i] = param2(i);
                bvec[i] = bvec2(i);
                dvec[i] = dvec2(i);
            }

            break;
        }
    }
    if (iter_lasso >= maxiter) {
        std::cout << "  Convergence NOT achieved within " << iter_lasso << " split-Bregman iterations." << std::endl;
    } else {
        std::cout << "  Convergence achieved in " << iter_lasso + 1 << " iterations." << std::endl;
    }

#else

    // Memory allocation

    double *Amat_1d;
    double *bvec2;
    double *dvec2;
    double *param_tmp;
    double *res;
    double **mat_precondition;
    double *vec_d;
    double *Qmat_CG;
    double *bvec_CG;
    double *bvec_CG_update;

    double dble_one = 1.0;
    double dble_zero = 0.0;
    double param_norm;

    unsigned long k;
    int N_ = N;
    int M_ = M;

    allocate(bvec2, N);
    allocate(dvec2, N);
    allocate(param_tmp, N);
    allocate(Amat_1d, M * N);
    allocate(Qmat_CG, N * N);
    allocate(bvec_CG, N);
    allocate(bvec_CG_update, N);
    allocate(res, M);

    // Prepare Qmat = A^T A, which is positive semidefinite.
    // When alpha is nonzero, Qmat becomes positive definite.

    k = 0;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j) {
            Amat_1d[k++] = Amat[j][i];
        }
    }
    dgemm_("T", "N", &N_, &N_, &M_, &dble_one, Amat_1d, &M_, Amat_1d, &M_, &dble_zero, Qmat_CG, &N_);
    for (i = 0; i < N; ++i) {
        Qmat_CG[(N + 1) * i] += alpha * alpha * lambda;
    }

    int ONE = 1;
    dgemv_("T", &M_, &N_, &dble_one, Amat_1d, &M_, fvec, &ONE, &dble_zero, bvec_CG, &ONE);

    if (preconditioner == 1) {
        //
        // Incomplete Cholesky factorization without fill-in, known as IC(0).
        //
        std::cout << "    LASSO_PCG = 1: Prepare preconditioning matrix ...";

        allocate(mat_precondition, N, N);
        allocate(vec_d, N);

        incomplete_cholesky_factorization(N, Qmat_CG, mat_precondition, vec_d);
        std::cout << " done." << std::endl << std::endl;

    }

    // Initialization

    if (cv_mode == 0) {
        std::cout << "    Start with b = 0 and d = 0" << std::endl;
        for (i = 0; i < N; ++i) {
            bvec2[i] = 0.0;
            dvec2[i] = 0.0;
            param_tmp[i] = 0.0;
        }
    } else if (cv_mode == 1) {
        std::cout << "    Start from b, d, and x of the previous run" << std::endl;
        for (i = 0; i < N; ++i) {
            bvec2[i] = bvec[i];
            dvec2[i] = dvec[i];
            param_tmp[i] = param[i];
        }
    } else if (cv_mode == 2) {
        std::cout << "    Start with b = 0 and d = 0" << std::endl;
        std::cout << "    Parameter x is initialized by solving the LS equation" << std::endl;
        for (i = 0; i < N; ++i) {
            bvec2[i] = 0.0;
            dvec2[i] = 0.0;
            param_tmp[i] = 0.0;
        }
        minimize_quadratic_CG(N, Qmat_CG, bvec_CG_update, param_tmp, 20 * N, true,
                              mat_precondition, vec_d, preconditioner);
    }


    for (iter_lasso = 0; iter_lasso < maxiter; ++iter_lasso) {

        do_print_log = !((iter_lasso + 1) % nfreq);

        if (do_print_log) {
            std::cout << "   SPLIT BREGMAN : " << std::setw(5) << iter_lasso + 1 << std::endl;
            std::cout << "    Update parameter vector using conjugate gradient" << std::endl;
        }

        // Update bvec_CG for CG

#pragma omp parallel for
        for (i = 0; i < N; ++i) {
            bvec_CG_update[i] = bvec_CG[i] + al * (dvec2[i] - bvec2[i]);
        }

        // Update parameters via CG minimization

        minimize_quadratic_CG(N, Qmat_CG, bvec_CG_update, param_tmp, maxiter_cg, do_print_log,
                              mat_precondition, vec_d, preconditioner);

        // Update dvec using the shrink operator
#pragma omp parallel for
        for (i = 0; i < N; ++i) {
            dvec2[i] = shrink(alpha * param_tmp[i] + bvec2[i], invlambda);
        }

        // Update bvec
#pragma omp parallel for
        for (i = 0; i < N; ++i) {
            bvec2[i] += alpha * param_tmp[i] - dvec2[i];
        }

        diff = 0.0;
#pragma omp parallel for reduction(+:diff)
        for (i = 0; i < N; ++i) {
            diff += std::pow(param[i] - param_tmp[i], 2);
        }

        if (do_print_log) {
            param_norm = 0.0;
            for (i = 0; i < N; ++i) {
                param_norm += std::pow(param_tmp[i], 2);
            }

            std::cout << std::endl;
            std::cout << "    1: ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << std::sqrt(
                    diff / static_cast<double>(N))
                << std::setw(15) << std::sqrt(diff / param_norm) << std::endl;

            tmp = 0.0;
            for (i = 0; i < N; ++i) {
                tmp += std::abs(param_tmp[i]);
            }
            std::cout << "    2: ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
            for (i = 0; i < M; ++i) {
                res[i] = -fvec[i];
            }
            dgemv_("N", &M_, &N_, &dble_one, Amat_1d, &M_, param_tmp, &ONE, &dble_one, res, &ONE);
            tmp = 0.0;
            for (i = 0; i < M; ++i) {
                tmp += res[i] * res[i];
            }
            std::cout << "    3: ||Au_{k}-f||_2          = " << std::setw(15)
                << std::sqrt(tmp) << std::setw(15) << std::sqrt(tmp / f2norm) << std::endl;
            tmp = 0.0;
            for (i = 0; i < N; ++i) {
                tmp += std::pow(dvec2[i] - alpha * param_tmp[i], 2.0);
            }
            std::cout << "    4: ||d_{k}-alpha*u_{k}||_2 = " << std::setw(15) << std::sqrt(tmp) << std::endl;
            tmp = 0.0;
            for (i = 0; i < N; ++i) {
                tmp += std::abs(dvec2[i]);
            }
            std::cout << "    5: ||d_{k}||_1             = " << std::setw(15) << tmp << std::endl;
            tmp = 0.0;
            for (i = 0; i < N; ++i) {
                tmp += bvec2[i] * bvec2[i];
            }
            std::cout << "    6: ||b_{k}||_2             = " << std::setw(15) << std::sqrt(tmp) << std::endl;
            std::cout << std::endl;

        }

        for (i = 0; i < N; ++i) {
            param[i] = param_tmp[i];
            bvec[i] = bvec2[i];
            dvec[i] = dvec2[i];
        }

        if (std::sqrt(diff / static_cast<double>(N)) < tolerance && iter_lasso >= miniter) {

            std::cout << "   SPLIT BREGMAN : " << std::setw(5) << iter_lasso + 1 << std::endl;
            std::cout << "    Update parameter vector using conjugate gradient" << std::endl;

            param_norm = 0.0;
            for (i = 0; i < N; ++i) {
                param_norm += std::pow(param_tmp[i], 2);
            }

            std::cout << std::endl;
            std::cout << "    1': ||u_{k}-u_{k-1}||_2     = "
                << std::setw(15) << std::sqrt(diff / static_cast<double>(N))
                << std::setw(15) << std::sqrt(diff / param_norm) << std::endl;

            tmp = 0.0;
            for (i = 0; i < N; ++i) {
                tmp += std::abs(param_tmp[i]);
            }
            std::cout << "    2': ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
            for (i = 0; i < M; ++i) {
                res[i] = -fvec[i];
            }
            dgemv_("N", &M_, &N_, &dble_one, Amat_1d, &M_, param_tmp, &ONE, &dble_one, res, &ONE);
            tmp = 0.0;
            for (i = 0; i < M; ++i) {
                tmp += res[i] * res[i];
            }
            std::cout << "    3': ||Au_{k}-f||_2          = " << std::setw(15)
                << std::sqrt(tmp) << std::setw(15) << std::sqrt(tmp / f2norm) << std::endl;
            tmp = 0.0;
            for (i = 0; i < N; ++i) {
                tmp += std::pow(dvec2[i] - alpha * param_tmp[i], 2.0);
            }
            std::cout << "    4': ||d_{k}-alpha*u_{k}||_2 = "
                << std::setw(15) << std::sqrt(tmp) << std::endl;
            tmp = 0.0;
            for (i = 0; i < N; ++i) {
                tmp += std::abs(dvec2[i]);
            }
            std::cout << "    5': ||d_{k}||_1             = "
                << std::setw(15) << tmp << std::endl;
            tmp = 0.0;
            for (i = 0; i < N; ++i) {
                tmp += bvec2[i] * bvec2[i];
            }
            std::cout << "    6': ||b_{k}||_2             = "
                << std::setw(15) << std::sqrt(tmp) << std::endl;
            std::cout << std::endl;

            for (i = 0; i < N; ++i) {
                param[i] = param_tmp[i];
                bvec[i] = bvec2[i];
                dvec[i] = dvec2[i];
            }

            break;
        }
    }
    if (iter_lasso >= maxiter) {
        std::cout << "  Convergence NOT achieved within "
            << iter_lasso << " split-Bregman iterations." << std::endl;
    } else {
        std::cout << "  Convergence achieved in " << iter_lasso + 1 << " iterations." << std::endl;
    }

    // Memory deallocation

    deallocate(bvec2);
    deallocate(dvec2);
    deallocate(param_tmp);
    deallocate(bvec_CG);
    deallocate(bvec_CG_update);
    deallocate(Qmat_CG);
    deallocate(Amat_1d);
    deallocate(res);
    if (preconditioner == 1) {
        deallocate(mat_precondition);
        deallocate(vec_d);
    }


#endif
}


void Lasso::minimize_quadratic_CG(const int N,
                                  double *Qmat,
                                  double *bvec,
                                  double *val_out,
                                  const int nmax_CG,
                                  const bool is_print_info,
                                  double **L,
                                  double *d,
                                  const int preconditioner)
{
    int i;
    int iter_CG;
    double tmp;
    double *descent_vec, *residual_vec;
    double *val_tmp, *vec_tmp;
    double alpha, beta;
    double diff, bnorm2;
    double dprod_res, dprod_res2;
    double dble_one = 1.0;
    double dble_minus = -1.0;
    double dble_zero = 0.0;
    int N_ = N;
    int ONE = 1;

    static double tol = 1.0e-10;

    allocate(descent_vec, N);
    allocate(val_tmp, N);
    allocate(residual_vec, N);
    allocate(vec_tmp, N);

    // Initialization

    bnorm2 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: bnorm2)
#endif
    for (i = 0; i < N; ++i) {
        val_tmp[i] = val_out[i];
        bnorm2 += bvec[i] * bvec[i];
        residual_vec[i] = bvec[i];
    }
    bnorm2 = 1.0 / bnorm2;

    dgemv_("N", &N_, &N_, &dble_minus, Qmat, &N_, val_tmp, &ONE, &dble_one, residual_vec, &ONE);

    if (preconditioner == 0) {

        for (i = 0; i < N; ++i) {
            descent_vec[i] = residual_vec[i];
        }

        for (iter_CG = 0; iter_CG < nmax_CG; ++iter_CG) {

            dgemv_("N", &N_, &N_, &dble_one, Qmat, &N_, descent_vec, &ONE, &dble_zero, vec_tmp, &ONE);

            dprod_res = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : dprod_res)
#endif
            for (i = 0; i < N; ++i) {
                dprod_res += residual_vec[i] * residual_vec[i];
            }
            tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tmp)
#endif
            for (i = 0; i < N; ++i) {
                tmp += descent_vec[i] * vec_tmp[i];
            }
            alpha = dprod_res / tmp;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
            for (i = 0; i < N; ++i) {
                val_tmp[i] += alpha * descent_vec[i];
                residual_vec[i] -= alpha * vec_tmp[i];
            }

            dprod_res2 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
            for (i = 0; i < N; ++i) {
                dprod_res2 += residual_vec[i] * residual_vec[i];
            }
            diff = dprod_res2 * bnorm2;

            if (is_print_info) {
                std::cout << "    CG " << std::setw(5) << iter_CG + 1 << ":";
                std::cout << " DIFF = " << std::sqrt(diff) << std::endl;
            }

            if (std::sqrt(diff) < tol) break;

            beta = dprod_res2 / dprod_res;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
            for (i = 0; i < N; ++i) {
                descent_vec[i] = residual_vec[i] + beta * descent_vec[i];
            }
        }
    } else if (preconditioner == 1) {

        forward_backward_substitution(N, L, d, residual_vec, descent_vec);
        dprod_res = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: dprod_res)
#endif
        for (i = 0; i < N; ++i) {
            dprod_res += residual_vec[i] * descent_vec[i];
        }

        for (iter_CG = 0; iter_CG < nmax_CG; ++iter_CG) {

            dgemv_("N", &N_, &N_, &dble_one, Qmat, &N_, descent_vec, &ONE, &dble_zero, vec_tmp, &ONE);

            tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: tmp)
#endif
            for (i = 0; i < N; ++i) {
                tmp += descent_vec[i] * vec_tmp[i];
            }
            alpha = dprod_res / tmp;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
            for (i = 0; i < N; ++i) {
                val_tmp[i] += alpha * descent_vec[i];
                residual_vec[i] -= alpha * vec_tmp[i];
            }
            diff = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: diff)
#endif
            for (i = 0; i < N; ++i) {
                diff += residual_vec[i] * residual_vec[i];
            }
            diff *= bnorm2;

            if (is_print_info) {
                std::cout << "    CG " << std::setw(5) << iter_CG + 1 << ":";
                std::cout << " DIFF = " << std::sqrt(diff) << std::endl;
            }
            if (std::sqrt(diff) < tol) break;

            forward_backward_substitution(N, L, d, residual_vec, vec_tmp);

            dprod_res2 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: dprod_res2)
#endif
            for (i = 0; i < N; ++i) {
                dprod_res2 += residual_vec[i] * vec_tmp[i];
            }
            beta = dprod_res2 / dprod_res;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (i = 0; i < N; ++i) {
                descent_vec[i] = vec_tmp[i] + beta * descent_vec[i];
            }
            dprod_res = dprod_res2;
        }

    } else {
        exit("minimize_quadratic_CG",
             "preconditioner > 1. This cannot happen.");
    }

    for (i = 0; i < N; ++i) {
        val_out[i] = val_tmp[i];
    }

    deallocate(descent_vec);
    deallocate(val_tmp);
    deallocate(residual_vec);
    deallocate(vec_tmp);
}


void Lasso::minimize_quadratic_CG(const int N,
                                  const Eigen::MatrixXd &Amat,
                                  const Eigen::VectorXd &bvec,
                                  Eigen::VectorXd &val_out,
                                  const int nmax_CG,
                                  const bool is_print_info,
                                  const Eigen::MatrixXd &L,
                                  const Eigen::VectorXd &d,
                                  const int preconditioner)
{
    int iter_CG;
    double alpha, beta;
    double diff;
    double bnorm2;
    double dprod_res, dprod_res2;
    static double tol = 1.0e-10;

    using namespace Eigen;

    VectorXd descent_vec(N);
    VectorXd vec_tmp(N), residual_vec(N);
    VectorXd val_tmp(N);

    // Initialization

    val_tmp = val_out;

    bnorm2 = 1.0 / bvec.dot(bvec);

    residual_vec = bvec - Amat * val_tmp;

    if (preconditioner == 0) {

        descent_vec = residual_vec;

        for (iter_CG = 0; iter_CG < nmax_CG; ++iter_CG) {

            vec_tmp = Amat * descent_vec;

            dprod_res = residual_vec.dot(residual_vec);
            alpha = dprod_res / descent_vec.dot(vec_tmp);

            val_tmp = val_tmp + alpha * descent_vec;
            residual_vec = residual_vec - alpha * vec_tmp;

            dprod_res2 = residual_vec.dot(residual_vec);

            diff = dprod_res2 * bnorm2;

            if (is_print_info) {
                std::cout << "    CG " << std::setw(5) << iter_CG + 1 << ":";
                std::cout << " DIFF = " << std::sqrt(diff) << std::endl;
            }

            if (std::sqrt(diff) < tol) break;

            beta = dprod_res2 / dprod_res;
            descent_vec = residual_vec + beta * descent_vec;
        }

    } else if (preconditioner == 1) {

        // p_0 = (LDL^T)^{-1} r_0
        forward_backward_substitution(N, L, d, residual_vec, descent_vec);
        // (r_0, p_0)
        dprod_res = residual_vec.dot(descent_vec);

        for (iter_CG = 0; iter_CG < nmax_CG; ++iter_CG) {

            // (A, p_i)
            vec_tmp = Amat * descent_vec;
            // Alpha = 
            alpha = dprod_res / descent_vec.dot(vec_tmp);

            val_tmp = val_tmp + alpha * descent_vec;
            residual_vec = residual_vec - alpha * vec_tmp;

            diff = residual_vec.dot(residual_vec) * bnorm2;
            //            diff = dprod_res2 * bnorm2;

            if (is_print_info) {
                std::cout << "    CG " << std::setw(5) << iter_CG + 1 << ":";
                std::cout << " DIFF = " << std::sqrt(diff) << std::endl;
            }
            if (std::sqrt(diff) < tol) break;

            forward_backward_substitution(N, L, d, residual_vec, vec_tmp);
            dprod_res2 = residual_vec.dot(vec_tmp);

            beta = dprod_res2 / dprod_res;
            descent_vec = vec_tmp + beta * descent_vec;

            dprod_res = dprod_res2;
        }

    } else {
        exit("minimize_quadratic_CG",
             "preconditioner > 1. This cannot happen.");
    }

    val_out = val_tmp;
}


void Lasso::calculate_residual(const int M,
                               const int N,
                               double **Amat,
                               double *param,
                               double *fvec,
                               const double fnorm,
                               double &res)
{
    int i, j;
    using namespace Eigen;

    MatrixXd Amat2(M, N);
    VectorXd param2(N);
    VectorXd fvec2(M), vec_tmp(M);

    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            Amat2(i, j) = Amat[i][j];
        }
        fvec2(i) = fvec[i];
    }

    for (i = 0; i < N; ++i) {
        param2(i) = param[i];
    }

    vec_tmp = Amat2 * param2 - fvec2;
    res = vec_tmp.dot(vec_tmp) / (fnorm * fnorm);
}

int Lasso::incomplete_cholesky_factorization(const int N,
                                             const Eigen::MatrixXd &A,
                                             Eigen::MatrixXd &L,
                                             Eigen::VectorXd &d)
{
    int i, j, k;
    double lld;
    double zero_criterion = 1.0e-8;

    if (N <= 0) return 0;

    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            L(i, j) = 0.0;
        }
    }

    L(0, 0) = A(0, 0);
    d(0) = 1.0 / L(0.0);

    for (i = 1; i < N; ++i) {
        for (j = 0; j <= i; ++j) {
            if (std::abs(A(i, j)) < zero_criterion) continue;

            lld = A(i, j);
            for (k = 0; k < j; ++k) {
                lld -= L(i, k) * L(j, k) * d(k);
            }
            L(i, j) = lld;
        }
        d(i) = 1.0 / L(i, i);
    }

    return 0;
}

int Lasso::incomplete_cholesky_factorization(const int N,
                                             double *A,
                                             double **L,
                                             double *d)
{
    // The vector A should be read as matrix A whose component is transposed (C<-->fortran convertion).

    double **mat;
    int i, j, k;
    unsigned long m;
    double lld;
    double zero_criterion = 1.0e-8;

    if (N <= 0) return 0;

    allocate(mat, N, N);

    m = 0;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            mat[j][i] = A[m++];
            L[i][j] = 0.0;
        }
    }

    L[0][0] = mat[0][0];
    d[0] = 1.0 / L[0][0];

    for (i = 1; i < N; ++i) {
        for (j = 0; j <= i; ++j) {
            if (std::abs(mat[i][j]) < zero_criterion) continue;

            lld = mat[i][j];
            for (k = 0; k < j; ++k) {
                lld -= L[i][k] * L[j][k] * d[k];
            }
            L[i][j] = lld;
        }
        d[i] = 1.0 / L[i][i];
    }

    deallocate(mat);

    return 0;
}


void Lasso::forward_backward_substitution(const int N,
                                          const Eigen::MatrixXd &L,
                                          const Eigen::VectorXd &d,
                                          const Eigen::VectorXd &vec_in,
                                          Eigen::VectorXd &vec_out)
{
    int i, j;
    double tmp;
    Eigen::VectorXd vec_tmp(N);

    for (i = 0; i < N; ++i) {
        vec_tmp(i) = 0.0;
        vec_out(i) = 0.0;
    }

    for (i = 0; i < N; ++i) {
        tmp = vec_in(i);
        for (j = 0; j < i; ++j) {
            tmp -= L(i, j) * vec_tmp(j);
        }
        vec_tmp(i) = tmp / L(i, i);
    }

    for (i = N - 1; i >= 0; --i) {
        tmp = 0.0;
        for (j = i + 1; j < N; ++j) {
            tmp += L(j, i) * vec_out(j);
        }
        vec_out(i) = vec_tmp(i) - d(i) * tmp;
    }
}


void Lasso::forward_backward_substitution(const int N,
                                          double **L,
                                          double *d,
                                          double *vec_in,
                                          double *vec_out)
{
    int i, j;
    double tmp;
    double *vec_tmp;

    allocate(vec_tmp, N);

    for (i = 0; i < N; ++i) {
        vec_tmp[i] = 0.0;
        vec_out[i] = 0.0;
    }

    for (i = 0; i < N; ++i) {
        tmp = vec_in[i];
        for (j = 0; j < i; ++j) {
            tmp -= L[i][j] * vec_tmp[j];
        }
        vec_tmp[i] = tmp / L[i][i];
    }

    for (i = N - 1; i >= 0; --i) {
        tmp = 0.0;
        for (j = i + 1; j < N; ++j) {
            tmp += L[j][i] * vec_out[j];
        }
        vec_out[i] = vec_tmp[i] - d[i] * tmp;
    }

    deallocate(vec_tmp);
}


void Lasso::coordinate_descent(const int M,
                               const int N,
                               const double alpha,
                               const double tolerance,
                               const int warm_start,
                               const int maxiter,
                               Eigen::VectorXd &x,
                               const Eigen::MatrixXd &A,
                               const Eigen::VectorXd &b,
                               const Eigen::VectorXd &C,
                               bool *has_prod,
                               Eigen::MatrixXd &Prod,
                               Eigen::VectorXd &grad,
                               const double f2norm,
                               const int nfreq,
                               Eigen::VectorXd scale_beta,
                               const int standardize)
{
    int i, j;
    int iloop;
    double diff;
    Eigen::VectorXd beta(N), delta(N);
    Eigen::VectorXd res(N);
    bool do_print_log;

    if (warm_start) {
        for (i = 0; i < N; ++i) beta(i) = x[i];
    } else {
        for (i = 0; i < N; ++i) beta(i) = 0.0;
        grad = C;
    }

    double Minv = 1.0 / static_cast<double>(M);

    iloop = 0;
    if (standardize) {
        while (iloop < maxiter) {
            do_print_log = !((iloop + 1) % nfreq);

            if (do_print_log) {
                std::cout << "   Coordinate Descent : " << std::setw(5) << iloop + 1 << std::endl;
            }
            delta = beta;
            for (i = 0; i < N; ++i) {
                beta(i) = shrink(Minv * grad(i) + beta(i), alpha);
                delta(i) -= beta(i);
                if (std::abs(delta(i)) > 0.0) {
                    if (!has_prod[i]) {
                        for (j = 0; j < N; ++j) {
                            Prod(j, i) = A.col(j).dot(A.col(i));
                        }
                        has_prod[i] = true;
                    }
                    grad = grad + Prod.col(i) * delta(i);
                }
            }
            ++iloop;
            diff = std::sqrt(delta.dot(delta) / static_cast<double>(N));

            if (diff < tolerance) break;

            if (do_print_log) {
                double param2norm = beta.dot(beta);
                std::cout << "    1: ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << diff
                    << std::setw(15) << diff * std::sqrt(static_cast<double>(N) / param2norm) << std::endl;
                double tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
                for (i = 0; i < N; ++i) {
                    tmp += std::abs(beta(i));
                }
                std::cout << "    2: ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
                res = A * beta - b;
                tmp = res.dot(res);
                std::cout << "    3: ||Au_{k}-f||_2          = " << std::setw(15) << std::sqrt(tmp)
                    << std::setw(15) << std::sqrt(tmp / f2norm) << std::endl;
                std::cout << std::endl;
            }
        }
    } else {
        // Non-standardized version. Needs additional operations
        while (iloop < maxiter) {
            do_print_log = !((iloop + 1) % nfreq);

            if (do_print_log) {
                std::cout << "   Coordinate Descent : " << std::setw(5) << iloop + 1 << std::endl;
            }
            delta = beta;
            for (i = 0; i < N; ++i) {
                beta(i) = shrink(Minv * grad(i) + beta(i) / scale_beta(i), alpha) * scale_beta(i);
                delta(i) -= beta(i);
                if (std::abs(delta(i)) > 0.0) {
                    if (!has_prod[i]) {
                        for (j = 0; j < N; ++j) {
                            Prod(j, i) = A.col(j).dot(A.col(i));
                        }
                        has_prod[i] = true;
                    }
                    grad = grad + Prod.col(i) * delta(i);
                }
            }
            ++iloop;
            diff = std::sqrt(delta.dot(delta) / static_cast<double>(N));

            if (diff < tolerance) break;

            if (do_print_log) {
                double param2norm = beta.dot(beta);
                std::cout << "    1: ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << diff
                    << std::setw(15) << diff * std::sqrt(static_cast<double>(N) / param2norm) << std::endl;
                double tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
                for (i = 0; i < N; ++i) {
                    tmp += std::abs(beta(i));
                }
                std::cout << "    2: ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
                res = A * beta - b;
                tmp = res.dot(res);
                std::cout << "    3: ||Au_{k}-f||_2          = " << std::setw(15) << std::sqrt(tmp)
                    << std::setw(15) << std::sqrt(tmp / f2norm) << std::endl;
                std::cout << std::endl;
            }
        }
    }


    if (iloop >= maxiter) {
        std::cout << "WARNING: Convergence NOT achieved within " << maxiter
            << " coordinate descent iterations." << std::endl;
    } else {
        std::cout << "  Convergence achieved in " << iloop << " iterations." << std::endl;
    }

    double param2norm = beta.dot(beta);
    if (std::abs(param2norm) < eps) {
        std::cout << "    1': ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << 0.0
            << std::setw(15) << 0.0 << std::endl;
    } else {
        std::cout << "    1': ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << diff
            << std::setw(15) << diff * std::sqrt(static_cast<double>(N) / param2norm) << std::endl;
    }

    double tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
    for (i = 0; i < N; ++i) {
        tmp += std::abs(beta(i));
    }
    std::cout << "    2': ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
    res = A * beta - b;
    tmp = res.dot(res);
    std::cout << "    3': ||Au_{k}-f||_2          = " << std::setw(15) << std::sqrt(tmp)
        << std::setw(15) << std::sqrt(tmp / f2norm) << std::endl;
    std::cout << std::endl;

    for (i = 0; i < N; ++i) x[i] = beta(i);
}


void Lasso::get_prefactor_force(const int maxorder,
                                Fcs *fcs,
                                Constraint *constraint,
                                Fitting *fitting,
                                std::vector<double> &prefactor)
{
    int i, j;
    int ishift2 = 0;
    int iparam2 = 0;
    int inew2, iold2;
    int iold2_dup;

    int *ind;

    allocate(ind, maxorder + 1);
    for (i = 0; i < maxorder; ++i) {
        for (const auto &it : constraint->index_bimap[i]) {
            inew2 = it.left + iparam2;
            iold2 = it.right;
            iold2_dup = 0;
            for (j = 0; j < iold2; ++j) {
                iold2_dup += fcs->nequiv[i][j];
            }

            for (j = 0; j < i + 2; ++j) {
                ind[j] = fcs->fc_table[i][iold2_dup].elems[j];
            }
            prefactor[inew2] = fitting->gamma(i + 2, ind);
        }

        ishift2 += fcs->nequiv[i].size();
        iparam2 += constraint->index_bimap[i].size();
    }
    deallocate(ind);
}
