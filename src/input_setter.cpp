/*
 input_setter.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <string>
#include "input_setter.h"
#include "alm_core.h"
#include "memory.h"
#include "files.h"
#include "interaction.h"
#include "system.h"
#include "symmetry.h"
#include "fitting.h"
#include "constraint.h"
#include "patterndisp.h"

using namespace ALM_NS;

InputSetter::InputSetter() {}

InputSetter::~InputSetter() {}

void InputSetter::set_general_vars(
    ALMCore *alm,
    const std::string prefix,
    const std::string mode,
    const std::string str_disp_basis,
    const std::string str_magmom,
    const int nat,
    const int nkd,
    const int nsym,
    const int is_printsymmetry,
    const int is_periodic[3],
    const bool trim_dispsign_for_evenfunc,
    const bool lspin,
    const bool print_hessian,
    const int noncollinear,
    const int trevsym,
    const std::string *kdname,
    const double * const *magmom,
    const double tolerance)
{
    int i, j;

    alm->files->job_title = prefix;
    alm->mode = mode;
    alm->system->nat = nat;
    alm->system->nkd = nkd;
    alm->system->str_magmom = str_magmom;
    alm->symmetry->nsym = nsym;
    alm->symmetry->is_printsymmetry = is_printsymmetry;
    alm->symmetry->tolerance = tolerance;
    alm->memory->allocate(alm->system->kdname, nkd);
    alm->memory->allocate(alm->system->magmom, nat, 3);

    for (i = 0; i < nkd; ++i) {
        alm->system->kdname[i] = kdname[i];
    }
    for (i = 0; i < 3; ++i) {
        alm->interaction->is_periodic[i] = is_periodic[i];
    }
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            alm->system->magmom[i][j] = magmom[i][j];
        }
    }
    alm->system->lspin = lspin;
    alm->system->noncollinear = noncollinear;
    alm->symmetry->trev_sym_mag = trevsym;
    alm->print_hessian = print_hessian;

    if (mode == "suggest") {
        alm->displace->disp_basis = str_disp_basis;
        alm->displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
    }
}

void InputSetter::set_cell_parameter(ALMCore *alm,
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

void InputSetter::set_interaction_vars(ALMCore *alm,
                                       const int maxorder,
                                       const int *nbody_include)
{
    int i;

    alm->interaction->maxorder = maxorder;
    alm->memory->allocate(alm->interaction->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i) {
        alm->interaction->nbody_include[i] = nbody_include[i];
    }
}

void InputSetter::set_cutoff_radii(ALMCore *alm,
                                   const int maxorder,
                                   const int nkd,
                                   const double * const * const * rcs)
{
    int i, j, k;

    alm->memory->allocate(alm->interaction->rcs, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                alm->interaction->rcs[i][j][k] = rcs[i][j][k];
            }
        }
    }
}

void InputSetter::set_fitting_vars(ALMCore *alm,
                                   const int ndata,
                                   const int nstart,
                                   const int nend,
                                   const int nskip,
                                   const int nboot,
                                   const std::string dfile,
                                   const std::string ffile,
                                   const int multiply_data,
                                   const int constraint_flag,
                                   const std::string rotation_axis,
                                   const std::string fc2_file,
                                   const std::string fc3_file,
                                   const bool fix_harmonic,
                                   const bool fix_cubic)
{
    alm->system->ndata = ndata;
    alm->system->nstart = nstart;
    alm->system->nend = nend;
    alm->system->nskip = nskip;

    alm->fitting->nboot = nboot;
    alm->files->file_disp = dfile;
    alm->files->file_force = ffile;
    alm->symmetry->multiply_data = multiply_data;
    alm->constraint->constraint_mode = constraint_flag;
    alm->constraint->rotation_axis = rotation_axis;
    alm->constraint->fc2_file = fc2_file;
    alm->constraint->fix_harmonic = fix_harmonic;
    alm->constraint->fc3_file = fc3_file;
    alm->constraint->fix_cubic = fix_cubic;
}

void InputSetter::set_atomic_positions(ALMCore *alm,
                                       const int nat,
                                       const int *kd,
                                       const double * const *xeq)
{
    int i, j;

    alm->memory->allocate(alm->system->xcoord, nat, 3);
    alm->memory->allocate(alm->system->kd, nat);

    for (i = 0; i < nat; ++i) {

        alm->system->kd[i] = kd[i];

        for (j = 0; j < 3; ++j) {
            alm->system->xcoord[i][j] = xeq[i][j];
        }
    }
}
