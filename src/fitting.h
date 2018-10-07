/*
 fitting.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>
#ifdef WITH_SPARSE_SOLVER
#include <Eigen/SparseCore>
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> SpMat;
#endif

#include "constraint.h"
#include "symmetry.h"
#include "fcs.h"
#include "timer.h"
#include "files.h"
#include <Eigen/Dense>


namespace ALM_NS
{
    class OptimizerControl
    {
    public:
        // General optimization options
        int optimizer;
        int use_sparse_solver;
        int maxnum_iteration;
        double tolerance_iteration;
        int output_frequency;

        // Options related to L1-regularized optimization
        int standardize;
        double displacement_scaling_factor;
        int debiase_after_l1opt;

        // cross-validation related variables
        int cross_validation_mode;
        int nset_cross_validation;
        double l1_alpha;
        double l1_alpha_min;
        double l1_alpha_max;
        int num_l1_alpha;
        int save_solution_path;

        OptimizerControl()
        {
            optimizer = 0;
            use_sparse_solver = 0;
            maxnum_iteration = 10000;
            tolerance_iteration = 1.0e-7;
            output_frequency = 1000;
            standardize = 1;
            displacement_scaling_factor = 1.0;
            debiase_after_l1opt = 0;
            cross_validation_mode = 0;
            nset_cross_validation = 1;
            l1_alpha = 0.0;
            l1_alpha_min = 1.0e-4;
            l1_alpha_max = 1.0;
            num_l1_alpha = 1;
            save_solution_path = 0;
        }
        ~OptimizerControl() = default;

        OptimizerControl(const OptimizerControl &obj) = default;
        OptimizerControl& operator=(const OptimizerControl &obj) = default;

    };

    class Fitting
    {
    public:
        Fitting();
        ~Fitting();

        int optimize_main(const Symmetry *symmetry,
                          const Constraint *constraint,
                          const Fcs *fcs,
                          const int maxorder,
                          const unsigned int nat,
                          const int verbosity,
                          const std::string file_disp,
                          const std::string file_force,
                          Timer *timer);

        int ndata, nstart, nend;
        int skip_s, skip_e;

        double *params;
        double **u_in;
        double **f_in;
        int use_sparseQR;


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


        int ndata_test, nstart_test, nend_test;
        std::string dfile_test, ffile_test;


        void lasso_main(const Symmetry *symmetry,
                        const Interaction *interaction,
                        const Fcs *fcs,
                        const Constraint *constraint,
                        const unsigned int nat,
                        const Files *files,
                        const int verbosity,
                        Fitting *fitting,
                        Timer *timer);

        void set_optimizer_control(const OptimizerControl &);
        OptimizerControl& get_optimizer_control();

    private:

        OptimizerControl optcontrol;
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

        double gamma(const int,
                     const int *) const;


        // Moved from lasso.h
        void coordinate_descent(int,
                                int,
                                double,
                                double,
                                int,
                                int,
                                Eigen::VectorXd &,
                                const Eigen::MatrixXd &,
                                const Eigen::VectorXd &,
                                const Eigen::VectorXd &,
                                bool *,
                                Eigen::MatrixXd &,
                                Eigen::VectorXd &,
                                double,
                                int,
                                Eigen::VectorXd,
                                int) const;

        void calculate_residual(int,
                                int,
                                double **,
                                double *,
                                double *,
                                double,
                                double &) const;


        void get_prefactor_force(const int,
                                 const Fcs *,
                                 const Constraint *,
                                 const Fitting *,
                                 std::vector<double> &) const;
    };

    inline double shrink(const double x,
                         const double a)
    {
        double xabs = std::abs(x);
        double sign = static_cast<double>((0.0 < x) - (x < 0.0));
        return sign * std::max<double>(xabs - a, 0.0);
    }

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
    void dgemm_(const char *TRANSA,
                const char *TRANSB,
                int *M,
                int *N,
                int *K,
                double *ALPHA,
                double *A,
                int *LDA,
                double *B,
                int *LDB,
                double *BETA,
                double *C,
                int *LDC);

    void dgemv_(const char *trans,
                int *M,
                int *N,
                double *alpha,
                double *a,
                int *lda,
                double *x,
                int *incx,
                double *beta,
                double *y,
                int *incy);
    }
}
