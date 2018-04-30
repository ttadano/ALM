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

InputSetter::InputSetter() {}

InputSetter::~InputSetter() {}

void InputSetter::deallocator(ALM *alm)
{
    if (alm->system->kdname) {
        deallocate(alm->system->kdname);
    }
    alm->system->kdname = nullptr;
    if (alm->system->xcoord) {
        deallocate(alm->system->xcoord);
    }
    alm->system->xcoord = nullptr;
    if (alm->system->kd) {
        deallocate(alm->system->kd);
    }
    alm->system->kd = nullptr;
    if (alm->system->magmom) {
        deallocate(alm->system->magmom);
    }
    alm->system->magmom = nullptr;
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
                                   const int nat,
                                   const int nkd,
                                   const int printsymmetry,
                                   const int is_periodic[3],
                                   const bool trim_dispsign_for_evenfunc,
                                   const bool lspin,
                                   const bool print_hessian,
                                   const int noncollinear,
                                   const int trevsym,
                                   const std::string *kdname,
                                   const double * const *magmom,
                                   const double tolerance,
                                   const double tolerance_constraint)
{
    int i, j;

    alm->files->job_title = prefix;
    alm->mode = mode;
    alm->set_verbose(verbosity);
    alm->system->nat = nat;
    alm->system->nkd = nkd;
    alm->system->str_magmom = str_magmom;
    alm->symmetry->printsymmetry = printsymmetry;
    alm->symmetry->tolerance = tolerance;
    allocate(alm->system->kdname, nkd);
    allocate(alm->system->magmom, nat, 3);

    for (i = 0; i < nkd; ++i) {
        alm->system->kdname[i] = kdname[i];
    }
    for (i = 0; i < 3; ++i) {
        alm->system->is_periodic[i] = is_periodic[i];
    }
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            alm->system->magmom[i][j] = magmom[i][j];
        }
    }
    alm->system->lspin = lspin;
    alm->system->noncollinear = noncollinear;
    alm->system->trev_sym_mag = trevsym;
    alm->files->print_hessian = print_hessian;
    alm->constraint->tolerance_constraint = tolerance_constraint;

    if (mode == "suggest") {
        alm->displace->disp_basis = str_disp_basis;
        alm->displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
    }
}

void InputSetter::set_cell_parameter(ALM *alm,
                                     const double a,
                                     const double lavec_tmp[3][3])
{
    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            alm->system->lavec[i][j] = a * lavec_tmp[i][j];
        }
    }
}

void InputSetter::set_interaction_vars(ALM *alm,
                                       const int maxorder,
                                       const int *nbody_include)
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
                                   const int nkd,
                                   const double * const * const *rcs)
{
    int i, j, k;

    allocate(alm->interaction->cutoff_radii, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                alm->interaction->cutoff_radii[i][j][k] = rcs[i][j][k];
            }
        }
    }
}

void InputSetter::set_fitting_vars(ALM *alm,
                                   const int ndata,
                                   const int nstart,
                                   const int nend,
                                   const std::string dfile,
                                   const std::string ffile,
                                   const int constraint_flag,
                                   const std::string rotation_axis,
                                   const std::string fc2_file,
                                   const std::string fc3_file,
                                   const bool fix_harmonic,
                                   const bool fix_cubic,
                                   const int flag_sparse)
{
    alm->fitting->ndata = ndata;
    alm->fitting->nstart = nstart;
    alm->fitting->nend = nend;
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

void InputSetter::set_atomic_positions(ALM *alm,
                                       const int nat,
                                       const int *kd,
                                       const double * const *xeq)
{
    int i, j;

    allocate(alm->system->xcoord, nat, 3);
    allocate(alm->system->kd, nat);

    for (i = 0; i < nat; ++i) {

        alm->system->kd[i] = kd[i];

        for (j = 0; j < 3; ++j) {
            alm->system->xcoord[i][j] = xeq[i][j];
        }
    }
}
