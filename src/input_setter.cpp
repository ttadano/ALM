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
#include "lasso.h"
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

    nbody_include = nullptr;
    cutoff_radii = nullptr;
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
    if (nbody_include) {
        deallocate(nbody_include);
    }
    if (cutoff_radii) {
        deallocate(cutoff_radii);
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

void InputSetter::set_interaction_vars(const int maxorder_in,
                                       const int *nbody_include_in)
{
    maxorder = maxorder_in;
    if (nbody_include) {
        deallocate(nbody_include);
    }
    allocate(nbody_include, maxorder);
    for (unsigned int i = 0; i < maxorder; i++) {
        nbody_include[i] = nbody_include_in[i];
    }
}

void InputSetter::set_cutoff_radii(const int maxorder_in,
                                   const unsigned int nkd_in,
                                   const double * const * const * cutoff_radii_in)
{
    if (cutoff_radii) {
        deallocate(cutoff_radii);
    }
    allocate(cutoff_radii, maxorder_in, nkd_in, nkd_in);
    for (unsigned int i = 0; i < maxorder_in; i++) {
        for (unsigned int j = 0; j < nkd_in; j++) {
            for (unsigned int k = 0; k < nkd_in; k++) {
                cutoff_radii[i][j][k] = cutoff_radii_in[i][j][k];
            }
        }
    }
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
    alm->symmetry->set_print_symmetry(printsymmetry);
    alm->symmetry->set_tolerance(tolerance);

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
    alm->constraint->set_tolerance_constraint(tolerance_constraint);

    if (mode == "suggest") {
        alm->displace->disp_basis = str_disp_basis;
        alm->displace->set_trim_dispsign_for_evenfunc(trim_dispsign_for_evenfunc);
    }
}

void InputSetter::define(ALM *alm)
{
    alm->define(maxorder,
                nkd,
                nbody_include,
                cutoff_radii);
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
    alm->fitting->set_ndata(ndata);
    alm->fitting->set_nstart(nstart);
    alm->fitting->set_nend(nend);
    alm->fitting->set_skip_s(skip_s);
    alm->fitting->set_skip_e(skip_e);
    alm->fitting->set_use_sparseQR(flag_sparse);

    alm->files->file_disp = dfile;
    alm->files->file_force = ffile;
    alm->constraint->set_constraint_mode(constraint_flag);
    alm->constraint->set_rotation_axis(rotation_axis);
    alm->constraint->set_fc_file(2, fc2_file);
    alm->constraint->set_fix_harmonic(fix_harmonic);
    alm->constraint->set_fc_file(3, fc3_file);
    alm->constraint->set_fix_cubic(fix_cubic);
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
    alm->lasso->disp_norm = lasso_dnorm;
    alm->lasso->l1_alpha = lasso_alpha;
    alm->lasso->l1_alpha_min = lasso_min_alpha;
    alm->lasso->l1_alpha_max = lasso_max_alpha;
    alm->lasso->num_l1_alpha = lasso_num_alpha;
    alm->lasso->l2_lambda = lasso_lambda;
    alm->lasso->lasso_tol = lasso_tol;
    alm->lasso->maxiter = lasso_maxiter;
    alm->lasso->maxiter_cg = lasso_maxiter_cg;
    alm->lasso->lasso_pcg = lasso_pcg;
    alm->lasso->lasso_cv = lasso_cv;
    alm->lasso->lasso_cvset = lasso_cvset;
    alm->lasso->save_solution_path = save_solution_path;
    alm->lasso->output_frequency = lasso_freq;
    alm->lasso->debias_ols = debias_ols;
    alm->lasso->lasso_zero_thr = lasso_zero_thr;
    alm->lasso->ndata_test = ndata_test;
    alm->lasso->nstart_test = nstart_test;
    alm->lasso->nend_test = nend_test;
    alm->lasso->dfile_test = dfile_test;
    alm->lasso->ffile_test = ffile_test;
    alm->lasso->standardize = standardize;
    alm->lasso->lasso_algo = lasso_algo;
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

void InputSetter::set_geometric_structure(ALM *alm)
{
    alm->set_cell(nat, lavec, xcoord, kd, kdname);
    alm->set_periodicity(is_periodic);
    alm->set_magnetic_params(nat, magmom, lspin, noncollinear, trevsym, str_magmom);
}
