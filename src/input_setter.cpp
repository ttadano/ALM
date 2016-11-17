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

void InputSetter::deallocate(ALMCore *alm_core)
{
    if (alm_core->system->kdname) {
        alm_core->memory->deallocate(alm_core->system->kdname);
    }
    alm_core->system->kdname = NULL;
    if (alm_core->system->xcoord) {
        alm_core->memory->deallocate(alm_core->system->xcoord);
    }
    alm_core->system->xcoord = NULL;
    if (alm_core->system->kd) {
        alm_core->memory->deallocate(alm_core->system->kd);
    }
    alm_core->system->kd = NULL;
    if (alm_core->system->magmom) {
        alm_core->memory->deallocate(alm_core->system->magmom);
    }
    alm_core->system->magmom = NULL;
    if (alm_core->interaction->nbody_include) {
        alm_core->memory->deallocate(alm_core->interaction->nbody_include);
    }
    alm_core->interaction->nbody_include = NULL;
    if (alm_core->interaction->rcs) {
        alm_core->memory->deallocate(alm_core->interaction->rcs);
    }
    alm_core->interaction->rcs = NULL;
}

void InputSetter::set_general_vars(
    ALMCore *alm_core,
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

    alm_core->files->job_title = prefix;
    alm_core->mode = mode;
    alm_core->system->nat = nat;
    alm_core->system->nkd = nkd;
    alm_core->system->str_magmom = str_magmom;
    alm_core->symmetry->nsym = nsym;
    alm_core->symmetry->is_printsymmetry = is_printsymmetry;
    alm_core->symmetry->tolerance = tolerance;
    alm_core->memory->allocate(alm_core->system->kdname, nkd);
    alm_core->memory->allocate(alm_core->system->magmom, nat, 3);

    for (i = 0; i < nkd; ++i) {
        alm_core->system->kdname[i] = kdname[i];
    }
    for (i = 0; i < 3; ++i) {
        alm_core->interaction->is_periodic[i] = is_periodic[i];
    }
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            alm_core->system->magmom[i][j] = magmom[i][j];
        }
    }
    alm_core->system->lspin = lspin;
    alm_core->system->noncollinear = noncollinear;
    alm_core->symmetry->trev_sym_mag = trevsym;
    alm_core->print_hessian = print_hessian;

    if (mode == "suggest") {
        alm_core->displace->disp_basis = str_disp_basis;
        alm_core->displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
    }
}

void InputSetter::set_cell_parameter(ALMCore *alm_core,
                                     const double a,
                                     const double lavec_tmp[3][3])
{
    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            alm_core->system->lavec[i][j] = a * lavec_tmp[i][j];
        }
    }
}

void InputSetter::set_interaction_vars(ALMCore *alm_core,
                                       const int maxorder,
                                       const int *nbody_include)
{
    int i;

    alm_core->interaction->maxorder = maxorder;
    alm_core->memory->allocate(alm_core->interaction->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i) {
        alm_core->interaction->nbody_include[i] = nbody_include[i];
    }
}

void InputSetter::set_cutoff_radii(ALMCore *alm_core,
                                   const int maxorder,
                                   const int nkd,
                                   const double * const * const * rcs)
{
    int i, j, k;

    alm_core->memory->allocate(alm_core->interaction->rcs, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                alm_core->interaction->rcs[i][j][k] = rcs[i][j][k];
            }
        }
    }
}

void InputSetter::set_fitting_vars(ALMCore *alm_core,
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
    alm_core->system->ndata = ndata;
    alm_core->system->nstart = nstart;
    alm_core->system->nend = nend;
    alm_core->system->nskip = nskip;

    alm_core->fitting->nboot = nboot;
    alm_core->files->file_disp = dfile;
    alm_core->files->file_force = ffile;
    alm_core->symmetry->multiply_data = multiply_data;
    alm_core->constraint->constraint_mode = constraint_flag;
    alm_core->constraint->rotation_axis = rotation_axis;
    alm_core->constraint->fc2_file = fc2_file;
    alm_core->constraint->fix_harmonic = fix_harmonic;
    alm_core->constraint->fc3_file = fc3_file;
    alm_core->constraint->fix_cubic = fix_cubic;
}

void InputSetter::set_atomic_positions(ALMCore *alm_core,
                                       const int nat,
                                       const int *kd,
                                       const double * const *xeq)
{
    int i, j;

    alm_core->memory->allocate(alm_core->system->xcoord, nat, 3);
    alm_core->memory->allocate(alm_core->system->kd, nat);

    for (i = 0; i < nat; ++i) {

        alm_core->system->kd[i] = kd[i];

        for (j = 0; j < 3; ++j) {
            alm_core->system->xcoord[i][j] = xeq[i][j];
        }
    }
}
