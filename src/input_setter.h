/*
 input_setter.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

//#include "alm.h"
#include "alm.h"
#include <string>

namespace ALM_NS
{
    class InputSetter
    {
    public:
        InputSetter();
        ~InputSetter();

        void deallocator(ALM *alm) const;
        void set_general_vars(ALM *alm,
                              std::string prefix,
                              std::string mode,
                              int verbosity,
                              std::string str_disp_basis,
                              std::string str_magmom,
                              int nat_in,
                              int nkd_in,
                              int is_printsymmetry,
                              const int is_periodic[3],
                              bool trim_dispsign_for_evenfunc,
                              bool lspin,
                              bool print_hessian,
                              int noncollinear,
                              int trevsym,
                              const std::string *kdname_in,
                              const double * const *magmom,
                              double tolerance,
                              double tolerance_constraint);
        void set_cell_parameter(const double a,
                                const double lavec_in[3][3]);
        void set_interaction_vars(ALM *alm,
                                  int maxorder,
                                  const int *nbody_include) const;
        void set_cutoff_radii(ALM *alm,
                              int maxorder,
                              int nkd_in,
                              const double * const * const *rcs) const;
        void set_fitting_vars(ALM *alm,
                              int ndata,
                              int nstart,
                              int nend,
                              int skip_s,
                              int skip_e,
                              std::string dfile,
                              std::string ffile,
                              int constraint_flag,
                              std::string rotation_axis,
                              std::string fc2_file,
                              std::string fc3_file,
                              bool fix_harmonic,
                              bool fix_cubic,
                              int flag_sparse) const;
        void set_lasso_vars(ALM *alm,
                            double lasso_alpha,
                            double lasso_min_alpha,
                            double lasso_max_alpha,
                            int lasso_num_alpha,
                            double lasso_tol,
                            int lasso_maxiter,
                            int lasso_freq,
                            int lasso_algo,
                            int standardize,
                            double lasso_dnorm,
                            double lasso_lambda,
                            int lasso_maxiter_cg,
                            int lasso_pcg,
                            int lasso_cv,
                            int lasso_cvset,
                            int save_solution_path,
                            int debias_ols,
                            double lasso_zero_thr,
                            int ndata_test,
                            int nstart_test,
                            int nend_test,
                            std::string dfile_test,
                            std::string ffile_test) const;
        void set_atomic_positions(const int nat_in,
                                  const int *kd_in,
                                  const double (*xcoord_in)[3]);
        void set_geometric_structure(ALM *alm) const;

    private:
        int nat, nkd;
        int *kd;
        double lavec[3][3];
        double (*xcoord)[3]; // fractional coordinate
        std::string *kdname;
        int is_periodic[3];

        bool lspin;
        double (*magmom)[3];
        int noncollinear;
        int trevsym;
        std::string str_magmom;
    };
}
