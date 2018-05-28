/*
 lasso.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm.h"
#include <string>
#include <Eigen/Dense>

namespace ALM_NS
{
    class Lasso
    {
    public:
        Lasso();
        ~Lasso();

        double disp_norm;
        double l1_alpha;
        double l1_alpha_min, l1_alpha_max;
        double l2_lambda;
        double lasso_tol;
        double lasso_zero_thr;
        int lasso_cv;
        int lasso_cvset;
        int maxiter;
        int maxiter_cg;
        int lasso_pcg;
        int output_frequency;
        int num_l1_alpha;
        int lasso_algo;
        int standardize;

        int ndata_test, nstart_test, nend_test;
        std::string dfile_test, ffile_test;

        void lasso_main(ALM *);

    private:

        void set_default_variables();

        void split_bregman_minimization(const int, const int, const double, const double,
                                        const double, const int, double **, double *, const double, double *,
                                        double *, double *, const int, const int);

        void coordinate_descent(const int, const int, const double, const double, const int, const int,
                                Eigen::VectorXd &, const Eigen::MatrixXd &, const Eigen::VectorXd &,
                                const Eigen::VectorXd &, bool *, Eigen::MatrixXd &, Eigen::VectorXd &,
                                const double, const int, Eigen::VectorXd, const int);

        void calculate_residual(const int, const int, double **, double *, double *,
                                const double, double &);

        void minimize_quadratic_CG(const int, double *, double *, double *, const int, const bool,
                                   double **, double *, const int);

        void minimize_quadratic_CG(const int, const Eigen::MatrixXd &, const Eigen::VectorXd &,
                                   Eigen::VectorXd &, const int, const bool,
                                   const Eigen::MatrixXd &, const Eigen::VectorXd &, const int);

        int incomplete_cholesky_factorization(const int, const Eigen::MatrixXd &,
                                              Eigen::MatrixXd &, Eigen::VectorXd &);

        int incomplete_cholesky_factorization(const int, double *, double **, double *);

        void forward_backward_substitution(const int, const Eigen::MatrixXd &, const Eigen::VectorXd &,
                                           const Eigen::VectorXd &, Eigen::VectorXd &);

        void forward_backward_substitution(const int, double **, double *, double *, double *);
    };

    inline double shrink(const double x, const double a)
    {
        double xabs = std::abs(x);
        double sign = static_cast<double>((0.0 < x) - (x < 0.0));
        return sign * std::max<double>(xabs - a, 0.0);
    }

    extern "C"
    {
        void dgemm_(const char *TRANSA, const char *TRANSB, int *M, int *N, int *K,
                    double *ALPHA, double *A, int *LDA, double *B, int *LDB, double *BETA,
                    double *C, int *LDC);

        void dgemv_(const char *trans, int *M, int *N, double *alpha, double *a, int *lda,
                    double *x, int *incx, double *beta, double *y, int *incy);
    }
}