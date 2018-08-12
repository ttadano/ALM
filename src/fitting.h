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

#ifdef WITH_SPARSE_SOLVER
#include <Eigen/SparseCore>
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> SpMat;
#endif

namespace ALM_NS
{
    class Fitting
    {
    public:
        Fitting();
        ~Fitting();

        int fitmain(ALM *);

        int ndata, nstart, nend;
        int skip_s, skip_e;

        double *params;
        double **u_in;
        double **f_in;
        int use_sparseQR;

        void set_displacement_and_force(const double * const *,
                                        const double * const *,
                                        int,
                                        int);


        void get_matrix_elements_algebraic_constraint(int,
                                                      int,
                                                      double *,
                                                      double *,
                                                      double &,
                                                      Symmetry *,
                                                      Fcs *,
                                                      Constraint *);

        void set_fcs_values(const int,
                            double *,
                            std::vector<int> *,
                            Constraint *);


        int get_ndata_used() const;
        double gamma(int,
                     const int *);

    private:

        int ndata_used;

        void set_default_variables();
        void deallocate_variables();

        void data_multiplier(double **,
                             std::vector<std::vector<double>> &,
                             int,
                             Symmetry *);

        int inprim_index(int,
                         Symmetry *);

        int fit_without_constraints(int,
                                    int,
                                    double *,
                                    double *,
                                    double *,
                                    int);

        int fit_algebraic_constraints(int,
                                      int,
                                      double *,
                                      double *,
                                      std::vector<double> &,
                                      double,
                                      int,
                                      Fcs *,
                                      Constraint *,
                                      int);

        int fit_with_constraints(int,
                                 int,
                                 int,
                                 double *,
                                 double *,
                                 double *,
                                 double **,
                                 double *,
                                 int);


        void get_matrix_elements(int,
                                 int,
                                 double *,
                                 double *,
                                 Symmetry *,
                                 Fcs *);

#ifdef WITH_SPARSE_SOLVER
        void get_matrix_elements_in_sparse_form(int,
                                                int,
                                                SpMat &,
                                                Eigen::VectorXd &,
                                                double &,
                                                Symmetry *,
                                                Fcs *,
                                                Constraint *);

        int run_eigen_sparseQR(const SpMat &,
                               const Eigen::VectorXd &,
                               std::vector<double> &, 
                               double,
                               int,
                               Fcs *,
                               Constraint *,
                               int);                          
#endif

        void recover_original_forceconstants(int,
                                             const std::vector<double> &,
                                             std::vector<double> &,
                                             std::vector<int> *,
                                             Constraint *);

        int factorial(int);
        int rankQRD(int,
                    int,
                    double *,
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
