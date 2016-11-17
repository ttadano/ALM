/*
 fitting.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <vector>
#include <set>
#include <string>
#ifdef _VSL
#include "mkl_vsl.h"
#endif

namespace ALM_NS
{
    class Fitting: protected Pointers
    {
    public:
        Fitting(class ALMCore *);
        ~Fitting();

        void fitmain();

        double *params;
        unsigned int seed;

        double **u_in;
        double **f_in;

        void set_displacement_and_force(const double * const * u_in,
                                        const double * const * f_in,
                                        const int nat,
                                        const int ndata_used);
        void calc_matrix_elements_algebraic_constraint(const int, const int, const int, const int,
                                                       const int, const int, const int, const int,
                                                       double **, double **, double **, double *, double *);

#ifdef _VSL
        VSLStreamStatePtr stream;
        int brng;
#endif

        double gamma(const int, const int *);

    private:

        void data_multiplier(double **u,
                             double **f,
                             const int nat,
                             const int ndata_used,
                             const int nmulti,
                             const int multiply_data);
        int get_number_for_multiplier(const int multiply_data);
        int inprim_index(const int);
        void fit_without_constraints(int, int, double **, double *, double *);
        void fit_algebraic_constraints(int, int, double **, double *,
                                       double *, double *, const int);

        void fit_with_constraints(int, int, int, double **, double *,
                                  double *, double **, double *);

        void calc_matrix_elements(const int, const int, const int,
                                  const int, const int, const int, const int,
                                  double **, double **, double **, double *);

        int factorial(const int);
        int rankSVD(const int, const int, double *, const double);
        int rankQRD(const int, const int, double *, const double);
        int rankSVD2(const int, const int, double **, const double);
#ifdef _USE_EIGEN_DISABLED
        int getRankEigen(const int, const int, double **);
#endif
    };

    extern "C"
    {
        void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda,
                     double *b, int *ldb, double *s, double *rcond, int *rank,
                     double *work, int *lwork, int *info);

        void dgglse_(int *m, int *n, int *p, double *a, int *lda,
                     double *b, int *ldb, double *c, double *d, double *x,
                     double *work, int *lwork, int *info);

        void dgesdd_(const char *jobz, int *m, int *n, double *a, int *lda,
                     double *s, double *u, int *ldu, double *vt, int *ldvt, double *work,
                     int *lwork, int *iwork, int *info);

        void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau,
                     double *work, int *lwork, int *info);

        void dgeqp3_(int *m, int *n, double *a, int *lda, int *jpvt,
                     double *tau, double *work, int *lwork, int *info);
    }
}

