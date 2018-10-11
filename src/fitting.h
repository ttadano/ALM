/*
 fitting.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>

#ifdef _VSL
#include "mkl_vsl.h"
#endif

#ifdef WITH_SPARSE_SOLVER
#include <Eigen/SparseCore>
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> SpMat;
#endif

#include "constraint.h"
#include "symmetry.h"
#include "fcs.h"
#include "timer.h"

namespace ALM_NS
{
    class Fitting
    {
    public:
        Fitting();
        ~Fitting();

        int fitmain(const Symmetry *symmetry,
                    const Constraint *constraint,
                    const Fcs *fcs,
                    const int maxorder,
                    const unsigned int nat,
                    const int verbosity,
                    const std::string file_disp,
                    const std::string file_force,
                    Timer *timer);

        void set_displacement_and_force(const double * const *,
                                        const double * const *,
                                        const int,
                                        const int);


        void get_matrix_elements_algebraic_constraint(const int,
                                                      const int,
                                                      double *,
                                                      double *,
                                                      double &,
                                                      const Symmetry *,
                                                      const Fcs *,
                                                      const Constraint *) const;

        void set_fcs_values(const int,
                            double *,
                            std::vector<int> *,
                            const Constraint *);


        int get_ndata_used() const;
        double gamma(const int,
                     const int *) const;

        int get_ndata() const;
        void set_ndata(const int);
        int get_nstart() const;
        void set_nstart(const int);
        int get_nend() const;
        void set_nend(const int);
        int get_skip_s() const;
        void set_skip_s(const int);
        int get_skip_e() const;
        void set_skip_e(const int);
        double *get_params() const;
        int get_use_sparseQR() const;
        void set_use_sparseQR(const int);

    private:

        int ndata, nstart, nend;
        int skip_s, skip_e;
        double *params;
        int use_sparseQR;
        double **u_in;
        double **f_in;

        int ndata_used;

        void set_default_variables();
        void deallocate_variables();

        void data_multiplier(double **,
                             std::vector<std::vector<double>> &,
                             const int,
                             const Symmetry *) const;

        int inprim_index(const int,
                         const Symmetry *) const;

        int fit_without_constraints(int,
                                    int,
                                    double *,
                                    double *,
                                    double *,
                                    const int) const;

        int fit_algebraic_constraints(int,
                                      int,
                                      double *,
                                      double *,
                                      std::vector<double> &,
                                      const double,
                                      const int,
                                      const Fcs *,
                                      const Constraint *,
                                      const int) const;

        int fit_with_constraints(int,
                                 int,
                                 int,
                                 double *,
                                 double *,
                                 double *,
                                 double **,
                                 double *,
                                 const int) const;


        void get_matrix_elements(const int,
                                 const int,
                                 double *,
                                 double *,
                                 const Symmetry *,
                                 const Fcs *) const;

#ifdef WITH_SPARSE_SOLVER
        void get_matrix_elements_in_sparse_form(const int,
                                                const int,
                                                SpMat &,
                                                Eigen::VectorXd &,
                                                double &,
                                                const Symmetry *,
                                                const Fcs *,
                                                const Constraint *);

        int run_eigen_sparseQR(const SpMat &,
                               const Eigen::VectorXd &,
                               std::vector<double> &,
                               const double,
                               const int,
                               const Fcs *,
                               const Constraint *,
                               const int);
#endif

        void recover_original_forceconstants(const int,
                                             const std::vector<double> &,
                                             std::vector<double> &,
                                             std::vector<int> *,
                                             const Constraint *) const;

        int factorial(const int) const;
        int rankQRD(const int,
                    const int,
                    double *,
                    const double) const;

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
