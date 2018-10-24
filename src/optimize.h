/*
 optimize.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>
#include "files.h"
#ifdef WITH_SPARSE_SOLVER
#include <Eigen/SparseCore>
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
#endif

#include "constraint.h"
#include "symmetry.h"
#include "fcs.h"
#include "timer.h"
#include <Eigen/Dense>


namespace ALM_NS
{
    class OptimizerControl
    {
    public:
        // General optimization options
        int optimizer; // 1 : least-squares, 2 : elastic net
        int use_sparse_solver; // 0: No, 1: Yes
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
        double l1_alpha; // L1-regularization coefficient
        double l1_alpha_min;
        double l1_alpha_max;
        int num_l1_alpha;
        double l1_ratio; // l1_ratio = 1 for LASSO; 0 < l1_ratio < 1 for Elastic net
        int save_solution_path;

        OptimizerControl()
        {
            optimizer = 1;
            use_sparse_solver = 0;
            maxnum_iteration = 10000;
            tolerance_iteration = 1.0e-15;
            output_frequency = 1000;
            standardize = 1;
            displacement_scaling_factor = 1.0;
            debiase_after_l1opt = 0;
            cross_validation_mode = 0;
            nset_cross_validation = 1;
            l1_alpha = 0.0;
            l1_alpha_min = 1.0e-4;
            l1_alpha_max = 1.0;
            l1_ratio = 1.0;
            num_l1_alpha = 1;
            save_solution_path = 0;
        }

        ~OptimizerControl() = default;

        OptimizerControl(const OptimizerControl &obj) = default;
        OptimizerControl& operator=(const OptimizerControl &obj) = default;
    };

    class Optimize
    {
    public:
        Optimize();
        ~Optimize();

        int optimize_main(const Symmetry *symmetry,
                          Constraint *constraint,
                          const Fcs *fcs,
                          const int maxorder,
                          const std::string file_prefix,
                          const std::vector<std::string> &str_order,
                          const unsigned int nat,
                          const int verbosity,
                          const DispForceFile &filedat_train,
                          const DispForceFile &filedat_test,
                          Timer *timer);

        //void set_displacement_and_force(const double * const *,
        //                                const double * const *,
        //                                const int,
        //                                const int);

        void set_training_data(const std::vector<std::vector<double>> &u_train_in,
                               const std::vector<std::vector<double>> &f_train_in);

        void set_test_data(const std::vector<std::vector<double>> &u_test_in,
                           const std::vector<std::vector<double>> &f_test_in);


        void get_matrix_elements_algebraic_constraint(const int,
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
        size_t get_number_of_rows_sensing_matrix() const;
        double* get_params() const;

        void set_optimizer_control(const OptimizerControl &);
        OptimizerControl get_optimizer_control() const;

        //int ndata_test, nstart_test, nend_test;
        //std::string dfile_test, ffile_test;

    private:

        double *params;
        /*    double **u_in;
            double **f_in;*/

        std::vector<std::vector<double>> u_train, f_train;
        std::vector<std::vector<double>> u_test, f_test;

        OptimizerControl optcontrol;
  //      int ndata_used;

        void set_default_variables();
        void deallocate_variables();

        void data_multiplier(const std::vector<std::vector<double>> &,
                             std::vector<std::vector<double>> &,
                             const Symmetry *) const;

        int inprim_index(const int,
                         const Symmetry *) const;

        int least_squares(const int maxorder,
                          const int natmin,
                          const int ntran,
                          const int N,
                          const int N_new,
                          const int M,
                          const int verbosity,
                          const Symmetry *symmetry,
                          const Fcs *fcs,
                          const Constraint *constraint,
                          std::vector<double> &param_out);

        int elastic_net(const std::string job_prefix,
                        const int maxorder,
                        const int natmin,
                        const int ntran,
                        const int N,
                        const int N_new,
                        const int M,
                        const int M_test,
                        const int ndata_used,
                        const int ndata_used_test,
                        const Symmetry *symmetry,
                        const std::vector<std::string> &str_order,
                        const Fcs *fcs,
                        Constraint *constraint,
                        const unsigned int nat,
                        const int verbosity,
                        std::vector<double> &param_out);


        int run_elastic_net_crossvalidation(const std::string job_prefix,
                                            const int maxorder,
                                            const int M,
                                            const int M_test,
                                            const int N_new,
                                            std::vector<double> &amat_1D,
                                            std::vector<double> &bvec,
                                            const double fnorm,
                                            std::vector<double> &amat_1D_test,
                                            std::vector<double> &bvec_test,
                                            const double fnorm_test,
                                            const Constraint *constraint,
                                            const int verbosity,
                                            std::vector<double> &param_out);

        int run_elastic_net_optimization(const int maxorder,
                                         const int M,
                                         const int N_new,
                                         std::vector<double> &amat_1D,
                                         std::vector<double> &bvec,
                                         const double fnorm,
                                         const std::vector<std::string> &str_order,
                                         const int verbosity,
                                         std::vector<double> &param_out);

        int run_least_squares_with_nonzero_coefs(const Eigen::MatrixXd &A_in,
                                                 const Eigen::VectorXd &b_in,
                                                 const Eigen::VectorXd &factor_std,
                                                 std::vector<double> &params,
                                                 const int verbosity);

        void get_standardizer(const Eigen::MatrixXd &Amat,
                              Eigen::VectorXd &mean,
                              Eigen::VectorXd &dev,
                              Eigen::VectorXd &factor_std,
                              Eigen::VectorXd &scale_beta);

        void apply_standardizer(Eigen::MatrixXd &Amat,
                                const Eigen::VectorXd &mean,
                                const Eigen::VectorXd &dev);

        double get_esimated_max_alpha(const Eigen::MatrixXd &Amat,
                                      const Eigen::VectorXd &bvec) const;


        int fit_without_constraints(const int,
                                    const int,
                                    double *,
                                    const double *,
                                    double *,
                                    const int) const;

        int fit_algebraic_constraints(const int,
                                      const int,
                                      double *,
                                      const double *,
                                      std::vector<double> &,
                                      const double,
                                      const int,
                                      const Fcs *,
                                      const Constraint *,
                                      const int) const;

        int fit_with_constraints(const int,
                                 const int,
                                 const int,
                                 double *,
                                 const double *,
                                 double *,
                                 const double * const *,
                                 double *,
                                 const int) const;


        void get_matrix_elements(const int,
                                 double *,
                                 double *,
                                 const Symmetry *,
                                 const Fcs *) const;

#ifdef WITH_SPARSE_SOLVER
        void get_matrix_elements_in_sparse_form(const int,
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
                                             const std::vector<int> *,
                                             const Constraint *) const;

        int factorial(const int) const;
        int rankQRD(const int,
                    const int,
                    double *,
                    const double) const;

        double gamma(const int,
                     const int *) const;


        // Moved from lasso.h
        void coordinate_descent(const int M,
                                const int N,
                                const double l1_alpha,
                                const int warm_start,
                                Eigen::VectorXd &x,
                                const Eigen::MatrixXd &A,
                                const Eigen::VectorXd &b,
                                const Eigen::VectorXd &grad0,
                                bool *has_prod,
                                Eigen::MatrixXd &Prod,
                                Eigen::VectorXd &grad,
                                const double fnorm,
                                const Eigen::VectorXd &scale_beta,
                                const int verbosity) const;
    };

    inline double shrink(const double x,
                         const double a)
    {
        const auto xabs = std::abs(x);
        const auto sign = static_cast<double>((0.0 < x) - (x < 0.0));
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
    }
}
