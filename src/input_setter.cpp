/*
 input_setter.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "input_setter.h"
// #include <iostream>
// #include <iomanip>
#include <string>
#include "memory.h"
#include "files.h"
#include "interaction.h"
#include "system.h"
#include "symmetry.h"
#include "error.h"
#include "fitting.h"
#include "constraint.h"
#include "patterndisp.h"

using namespace ALM_NS;

InputSetter::InputSetter(ALMCore *alm): Pointers(alm) {}

InputSetter::~InputSetter() {}

void InputSetter::set_input(const std::string prefix,
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
                            const double tolerance,
                            const double a,
                            const double lavec_tmp[3][3],
                            const int maxorder,
                            const int *nbody_include,
                            const double * const * const * rcs,
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
                            const bool fix_cubic,
                            const int *kd,
                            const double * const *xeq)
{
    set_general_vars(prefix,
                     mode,
                     str_disp_basis,
                     str_magmom,
                     nat,
                     nkd,
                     nsym,
                     is_printsymmetry,
                     is_periodic,
                     trim_dispsign_for_evenfunc,
                     lspin,
                     print_hessian,
                     noncollinear,
                     trevsym,
                     kdname,
                     magmom,
                     tolerance);
    set_cell_parameter(a, lavec_tmp);
    set_interaction_vars(maxorder, nbody_include);
    set_cutoff_radii(maxorder, nkd, rcs);
    set_fitting_vars(ndata,
                     nstart,
                     nend,
                     nskip,
                     nboot,
                     dfile,
                     ffile,
                     multiply_data,
                     constraint_flag,
                     rotation_axis,
                     fc2_file,
                     fc3_file,
                     fix_harmonic,
                     fix_cubic);
    set_atomic_positions(nat, kd, xeq);
}

void InputSetter::set_general_vars(
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

    files->job_title = prefix;
    alm->mode = mode;
    system->nat = nat;
    system->nkd = nkd;
    system->str_magmom = str_magmom;
    symmetry->nsym = nsym;
    symmetry->is_printsymmetry = is_printsymmetry;
    symmetry->tolerance = tolerance;
    memory->allocate(system->kdname, nkd);
    memory->allocate(system->magmom, nat, 3);

    for (i = 0; i < nkd; ++i) {
        system->kdname[i] = kdname[i];
    }
    for (i = 0; i < 3; ++i) {
        interaction->is_periodic[i] = is_periodic[i];
    }
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            system->magmom[i][j] = magmom[i][j];
        }
    }
    system->lspin = lspin;
    system->noncollinear = noncollinear;
    symmetry->trev_sym_mag = trevsym;
    alm->print_hessian = print_hessian;

    if (mode == "suggest") {
        displace->disp_basis = str_disp_basis;
        displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
    }
}

void InputSetter::set_cell_parameter(const double a,
                                     const double lavec_tmp[3][3])
{
    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            system->lavec[i][j] = a * lavec_tmp[i][j];
        }
    }
}

void InputSetter::set_interaction_vars(const int maxorder,
                                       const int *nbody_include)
{
    int i;

    interaction->maxorder = maxorder;
    memory->allocate(interaction->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i) {
        interaction->nbody_include[i] = nbody_include[i];
    }
}

void InputSetter::set_cutoff_radii(const int maxorder,
                                   const int nkd,
                                   const double * const * const * rcs)
{
    int i, j, k;

    memory->allocate(interaction->rcs, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                interaction->rcs[i][j][k] = rcs[i][j][k];
            }
        }
    }
}

void InputSetter::set_fitting_vars(const int ndata,
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
    system->ndata = ndata;
    system->nstart = nstart;
    system->nend = nend;
    system->nskip = nskip;

    fitting->nboot = nboot;
    files->file_disp = dfile;
    files->file_force = ffile;
    symmetry->multiply_data = multiply_data;
    constraint->constraint_mode = constraint_flag;
    constraint->rotation_axis = rotation_axis;
    constraint->fc2_file = fc2_file;
    constraint->fix_harmonic = fix_harmonic;
    constraint->fc3_file = fc3_file;
    constraint->fix_cubic = fix_cubic;
}

void InputSetter::set_atomic_positions(const int nat,
                                       const int *kd,
                                       const double * const *xeq)
{
    int i, j;

    memory->allocate(system->xcoord, nat, 3);
    memory->allocate(system->kd, nat);

    for (i = 0; i < nat; ++i) {

        system->kd[i] = kd[i];

        for (j = 0; j < 3; ++j) {
            system->xcoord[i][j] = xeq[i][j];
        }
    }
}
