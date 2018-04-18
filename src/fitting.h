/*
 fitting.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm.h"
#include <vector>
#ifdef _VSL
#include "mkl_vsl.h"
#endif

namespace ALM_NS
{
    class Fitting
    {
    public:
        Fitting();
        ~Fitting();

        void fitmain(ALM *);

        int ndata, nstart, nend;

        double *params;
        double **u_in;
        double **f_in;

        void set_displacement_and_force(const double * const *u_in,
                                        const double * const *f_in,
                                        int nat,
                                        int ndata_used);
        void calc_matrix_elements_algebraic_constraint(int,
                                                       int,
                                                       int,
                                                       int,
                                                       int,
                                                       int,
                                                       int,
                                                       int,
                                                       double **,
                                                       double **,
                                                       double *,
                                                       double *,
                                                       double *,
                                                       Symmetry *,
                                                       Fcs *,
                                                       Constraint *);

        void get_matrix_elements(int,
                                 int,
                                 int,
                                 double *,
                                 double *,
                                 Symmetry *,
                                 Fcs *);

        const int get_ndata_used();
        double gamma(int, const int *);

    private:

        int ndata_used;

        void set_default_variables();
        void deallocate_variables();
        void data_multiplier(double **u,
                             double **f,
                             int nat,
                             int ndata_used,
                             Symmetry *);

        void data_multiplier(const std::vector<std::vector<double>> &,
                             std::vector<std::vector<double>> &,
                             int,
                             Symmetry *);

        int inprim_index(int, Symmetry *);

        void fit_without_constraints(int,
                                     int,
                                     double *,
                                     double *,
                                     double *);

        void fit_algebraic_constraints(int,
                                       int,
                                       double *,
                                       double *,
                                       std::vector<double> &,
                                       double,
                                       int,
                                       Fcs *,
                                       Constraint *);

        void fit_with_constraints(int,
                                  int,
                                  int,
                                  double *,
                                  double *,
                                  double *,
                                  double **,
                                  double *);

        void calc_matrix_elements(int,
                                  int,
                                  int,
                                  int,
                                  int,
                                  int,
                                  double **,
                                  double **,
                                  double *,
                                  double *,
                                  Symmetry *,
                                  Fcs *);


        void get_matrix_elements(int,
                                 int,
                                 const std::vector<std::vector<double>> &,
                                 const std::vector<std::vector<double>> &,
                                 std::vector<double> &,
                                 std::vector<double> &,
                                 Symmetry *,
                                 Fcs *);


        void get_matrix_elements_algebraic_constraint(int,
                                                      int,
                                                      const std::vector<std::vector<double>> &,
                                                      const std::vector<std::vector<double>> &,
                                                      std::vector<double> &,
                                                      std::vector<double> &,
                                                      double &,
                                                      Symmetry *,
                                                      Fcs *,
                                                      Constraint *);

        void recover_original_forceconstants(int,
                                             const std::vector<double> &,
                                             std::vector<double> &,
                                             std::vector<int> *,
                                             Constraint *);

        int factorial(int);
        int rankSVD(int,
                    int,
                    double *,
                    double);
        int rankQRD(int,
                    int,
                    double *,
                    double);
        int rankSVD2(int,
                     int,
                     double **,
                     double);
    };

    extern "C" {
    void dgelss_(int *m,
                 int *n,
                 int *nrhs,
                 double *a,
                 int *lda,
                 double *b,
                 int *ldb,
                 double *s,
                 double *rcond,
                 int *rank,
                 double *work,
                 int *lwork,
                 int *info);

    void dgglse_(int *m,
                 int *n,
                 int *p,
                 double *a,
                 int *lda,
                 double *b,
                 int *ldb,
                 double *c,
                 double *d,
                 double *x,
                 double *work,
                 int *lwork,
                 int *info);

    void dgesdd_(const char *jobz,
                 int *m,
                 int *n,
                 double *a,
                 int *lda,
                 double *s,
                 double *u,
                 int *ldu,
                 double *vt,
                 int *ldvt,
                 double *work,
                 int *lwork,
                 int *iwork,
                 int *info);

    void dgeqrf_(int *m,
                 int *n,
                 double *a,
                 int *lda,
                 double *tau,
                 double *work,
                 int *lwork,
                 int *info);

    void dgeqp3_(int *m,
                 int *n,
                 double *a,
                 int *lda,
                 int *jpvt,
                 double *tau,
                 double *work,
                 int *lwork,
                 int *info);
    }
}
