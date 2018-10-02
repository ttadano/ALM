/*
 input_setter.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <string>
#include "input_setter.h"
#include "memory.h"
#include "files.h"
#include "interaction.h"
#include "system.h"
#include "symmetry.h"
#include "fitting.h"
#include "constraint.h"
#include "patterndisp.h"
#include "alm.h"

using namespace ALM_NS;

InputSetter::InputSetter() {
    nat = 0;
    nkd = 0;
    kd = nullptr;
    kdname = nullptr;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; j++) {
            lavec[i][j] = 0.0;
        }
    }
    xcoord = nullptr;
    is_periodic[0] = 1;
    is_periodic[1] = 1;
    is_periodic[2] = 1;

    lspin = false;
    magmom = nullptr;
    noncollinear = 0;
    trevsym = 1;
    str_magmom = "";
}

InputSetter::~InputSetter() {
    if (kdname) {
        deallocate(kdname);
    }
    if (xcoord) {
        deallocate(xcoord);
    }
    if (kd) {
        deallocate(kd);
    }
    if (magmom) {
        deallocate(magmom);
    }
}

void InputSetter::deallocator(ALM *alm) const
{
    if (alm->interaction->nbody_include) {
        deallocate(alm->interaction->nbody_include);
    }
    alm->interaction->nbody_include = nullptr;
    if (alm->interaction->cutoff_radii) {
        deallocate(alm->interaction->cutoff_radii);
    }
    alm->interaction->cutoff_radii = nullptr;
}

void InputSetter::set_general_vars(ALM *alm,
                                   const std::string prefix,
                                   const std::string mode,
                                   const int verbosity,
                                   const std::string str_disp_basis,
                                   const std::string str_magmom,
                                   const int nat_in,
                                   const int nkd_in,
                                   const int printsymmetry,
                                   const int is_periodic_in[3],
                                   const bool trim_dispsign_for_evenfunc,
                                   const bool lspin_in,
                                   const bool print_hessian,
                                   const int noncollinear_in,
                                   const int trevsym_in,
                                   const std::string *kdname_in,
                                   const double * const *magmom_in,
                                   const double tolerance,
                                   const double tolerance_constraint)
{
    int i, j;

    alm->files->job_title = prefix;
    alm->set_run_mode(mode);
    alm->set_verbosity(verbosity);
    nat = nat_in;
    nkd = nkd_in;
    alm->symmetry->printsymmetry = printsymmetry;
    alm->symmetry->tolerance = tolerance;

    if (kdname) {
        deallocate(kdname);
    }
    allocate(kdname, nkd);
    for (i = 0; i < nkd; ++i) {
        kdname[i] = kdname_in[i];
    }

    if (magmom) {
        deallocate(magmom);
    }
    allocate(magmom, nat);

    for (i = 0; i < nat; i ++) {
        for (j = 0; j < 3; j++) {
            magmom[i][j] = magmom_in[i][j];
        }
    }
    lspin = lspin_in;
    noncollinear = noncollinear_in;
    trevsym = trevsym_in;

    for (i = 0; i < 3; i++) {
        is_periodic[i] = is_periodic_in[i];
    }

    alm->files->print_hessian = print_hessian;
    alm->constraint->tolerance_constraint = tolerance_constraint;

    if (mode == "suggest") {
        alm->displace->disp_basis = str_disp_basis;
        alm->displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
    }
}

void InputSetter::set_cell_parameter(const double a,
                                     const double lavec_in[3][3])
{
    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            lavec[i][j] = a * lavec_in[i][j];
        }
    }
}

void InputSetter::set_interaction_vars(ALM *alm,
                                       const int maxorder,
                                       const int *nbody_include) const
{
    int i;

    alm->interaction->maxorder = maxorder;
    allocate(alm->interaction->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i) {
        alm->interaction->nbody_include[i] = nbody_include[i];
    }
}

void InputSetter::set_cutoff_radii(ALM *alm,
                                   const int maxorder,
                                   const int nkd_in,
                                   const double * const * const *rcs) const
{
    int i, j, k;

    allocate(alm->interaction->cutoff_radii, maxorder, nkd_in, nkd_in);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd_in; ++j) {
            for (k = 0; k < nkd_in; ++k) {
                alm->interaction->cutoff_radii[i][j][k] = rcs[i][j][k];
            }
        }
    }
}

void InputSetter::set_fitting_vars(ALM *alm,
                                   const int ndata,
                                   const int nstart,
                                   const int nend,
                                   const int skip_s,
                                   const int skip_e,
                                   const std::string dfile,
                                   const std::string ffile,
                                   const int constraint_flag,
                                   const std::string rotation_axis,
                                   const std::string fc2_file,
                                   const std::string fc3_file,
                                   const bool fix_harmonic,
                                   const bool fix_cubic,
                                   const int flag_sparse) const
{
    alm->fitting->ndata = ndata;
    alm->fitting->nstart = nstart;
    alm->fitting->nend = nend;
    alm->fitting->skip_s = skip_s;
    alm->fitting->skip_e = skip_e;
    alm->fitting->use_sparseQR = flag_sparse;

    alm->files->file_disp = dfile;
    alm->files->file_force = ffile;
    alm->constraint->constraint_mode = constraint_flag;
    alm->constraint->rotation_axis = rotation_axis;
    alm->constraint->fc2_file = fc2_file;
    alm->constraint->fix_harmonic = fix_harmonic;
    alm->constraint->fc3_file = fc3_file;
    alm->constraint->fix_cubic = fix_cubic;
}

void InputSetter::set_lasso_vars(ALM *alm,
                                 const double lasso_alpha,
                                 const double lasso_min_alpha,
                                 const double lasso_max_alpha,
                                 const int lasso_num_alpha,
                                 const double lasso_tol,
                                 const int lasso_maxiter,
                                 const int lasso_freq,
                                 const int lasso_algo,
                                 const int standardize,
                                 const double lasso_dnorm,
                                 const double lasso_lambda,
                                 const int lasso_maxiter_cg,
                                 const int lasso_pcg,
                                 const int lasso_cv,
                                 const int lasso_cvset,
                                 const int save_solution_path,
                                 const int debias_ols,
                                 const double lasso_zero_thr,
                                 const int ndata_test,
                                 const int nstart_test,
                                 const int nend_test,
                                 std::string dfile_test,
                                 std::string ffile_test) const
{
    alm->fitting->disp_norm = lasso_dnorm;
    alm->fitting->l1_alpha = lasso_alpha;
    alm->fitting->l1_alpha_min = lasso_min_alpha;
    alm->fitting->l1_alpha_max = lasso_max_alpha;
    alm->fitting->num_l1_alpha = lasso_num_alpha;
    alm->fitting->l2_lambda = lasso_lambda;
    alm->fitting->lasso_tol = lasso_tol;
    alm->fitting->maxiter = lasso_maxiter;
    alm->fitting->maxiter_cg = lasso_maxiter_cg;
    alm->fitting->lasso_pcg = lasso_pcg;
    alm->fitting->lasso_cv = lasso_cv;
    alm->fitting->lasso_cvset = lasso_cvset;
    alm->fitting->save_solution_path = save_solution_path;
    alm->fitting->output_frequency = lasso_freq;
    alm->fitting->debias_ols = debias_ols;
    alm->fitting->lasso_zero_thr = lasso_zero_thr;
    alm->fitting->ndata_test = ndata_test;
    alm->fitting->nstart_test = nstart_test;
    alm->fitting->nend_test = nend_test;
    alm->fitting->dfile_test = dfile_test;
    alm->fitting->ffile_test = ffile_test;
    alm->fitting->standardize = standardize;
    alm->fitting->lasso_algo = lasso_algo;
}

void InputSetter::set_atomic_positions(const int nat_in,
                                       const int *kd_in,
                                       const double (*xcoord_in)[3])
{
    int i, j;

    if (kd) {
        deallocate(kd);
    }
    if (xcoord) {
        deallocate(xcoord);
    }
    allocate(xcoord, nat);
    allocate(kd, nat);

    for (i = 0; i < nat; ++i) {
        kd[i] = kd_in[i];
        for (j = 0; j < 3; ++j) {
            xcoord[i][j] = xcoord_in[i][j];
        }
    }
}

void InputSetter::set_geometric_structure(ALM *alm) const
{
    alm->set_cell(nat, lavec, xcoord, kd, kdname);
    alm->set_periodicity(is_periodic);
    alm->set_magnetic_params(magmom, lspin, noncollinear, trevsym, str_magmom);
}
