/*
 alm.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <string>
#include "alm.h"
#include "alm_core.h"
#include "constraint.h"
#include "fcs.h"
#include "fitting.h"
#include "files.h"
#include "interaction.h"
#include "memory.h"
#include "symmetry.h"
#include "system.h"
#include "patterndisp.h"

using namespace ALM_NS;

ALM::ALM()
{
    alm_core = new ALMCore();
    alm_core->create();
}

ALM::~ALM()
{
    delete alm_core;
}

void ALM::initialize()
{
    alm_core->initialize();
}

void ALM::finalize()
{
    alm_core->finalize();
}

void ALM::set_run_mode(const std::string mode)
{
    alm_core->mode = mode;
}

void ALM::set_output_control_params(const std::string prefix, // PREFIX
                                    const int is_printsymmetry, // PRINTSYM
                                    const bool print_hessian) // HESSIAN
{
    alm_core->files->job_title = prefix;
    alm_core->symmetry->is_printsymmetry = is_printsymmetry;
    alm_core->print_hessian = print_hessian;
}

void ALM::set_symmetry_params(const int nsym, // NSYM
                              const double tolerance) // TOLERANCE
{
    alm_core->symmetry->nsym = nsym;
    alm_core->symmetry->tolerance = tolerance;
}

void ALM::set_displacement_params(const std::string str_disp_basis, // DBASIS
                                  const bool trim_dispsign_for_evenfunc) // TRIMEVEN
{
    alm_core->displace->disp_basis = str_disp_basis;
    alm_core->displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
}

void ALM::set_periodicity(const int is_periodic[3]) // PERIODIC
{
    int i;
    for (i = 0; i < 3; ++i) {
        alm_core->interaction->is_periodic[i] = is_periodic[i];
    }
}

void ALM::set_cell(const int nat,
                   const int nkd,
                   const double lavec[3][3],
                   const double * const *xcoord,
                   const int *kd,
                   const std::string *kdname)
{
    int i, j;

    alm_core->system->nat = nat;
    alm_core->system->nkd = nkd;
    alm_core->memory->allocate(alm_core->system->kdname, nkd);
    alm_core->memory->allocate(alm_core->system->xcoord, nat, 3);
    alm_core->memory->allocate(alm_core->system->kd, nat);
    for (i = 0; i < nkd; ++i) {
        alm_core->system->kdname[i] = kdname[i];
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            alm_core->system->lavec[i][j] = lavec[i][j];
        }
    }
    for (i = 0; i < nat; ++i) {
        alm_core->system->kd[i] = kd[i];
        for (j = 0; j < 3; ++j) {
            alm_core->system->xcoord[i][j] = xcoord[i][j];
        }
    }
}

void ALM::set_magnetic_params(const double * const *magmom, // MAGMOM
                              const bool lspin, // MAGMOM
                              const int noncollinear, // NONCOLLINEAR
                              const int trevsym, // TREVSYM
                              const std::string str_magmom) // MAGMOM
{
    int i, j, nat;

    nat = alm_core->system->nat;
    alm_core->system->lspin = lspin;
    alm_core->system->noncollinear = noncollinear;
    alm_core->system->str_magmom = str_magmom;
    alm_core->memory->allocate(alm_core->system->magmom, nat, 3);
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            alm_core->system->magmom[i][j] = magmom[i][j];
        }
    }
}

void ALM::set_force_file_params(const int ndata, // NDATA
                                const int nstart, // NSTART
                                const int nend) // NEND
{
    alm_core->system->ndata = ndata;
    alm_core->system->nstart = nstart;
    alm_core->system->nend = nend;
}

void ALM::set_fitting_constraint(const int constraint_flag, // ICONST
                                 const std::string rotation_axis) // ROTAXIS
{
    alm_core->constraint->constraint_mode = constraint_flag;
    alm_core->constraint->rotation_axis = rotation_axis;
}

void ALM::set_fitting_params(const int nskip, // NSKIP
                             const int nboot, // NBOOT
                             const int multiply_data) // MULTDAT
{
    alm_core->system->nskip = nskip;
    alm_core->fitting->nboot = nboot;
    alm_core->symmetry->multiply_data = multiply_data;
}

void ALM::set_fitting_filenames(const std::string dfile, // DFILE
                                const std::string ffile) // FFILE
{
    alm_core->files->file_disp = dfile;
    alm_core->files->file_force = ffile;
}

void ALM::set_fitting_fc2_filename(const std::string fc2_file) // FC2XML
{
    bool fix_harmonic;

    alm_core->constraint->fc2_file = fc2_file;
    if (fc2_file.empty()) {
        fix_harmonic = false;
    } else {
        fix_harmonic = true;
    }
    alm_core->constraint->fix_harmonic = fix_harmonic;
}

void ALM::set_fitting_fc3_filename(const std::string fc3_file) // FC3XML
{
    bool fix_cubic;

    alm_core->constraint->fc3_file = fc3_file;
    if (fc3_file.empty()) {
        fix_cubic = false;
    } else {
        fix_cubic = true;
    }
    alm_core->constraint->fix_cubic = fix_cubic;
}

void ALM::set_interaction_vars(const int maxorder, // NORDER harmonic=1
                               const int *nbody_include) // NBODY
{
    int i;

    // nbody_include is defined as follows:
    //
    // if (nbody_v[0].empty()) { // Default [2, 3, 4, ..., NORDER + 1]
    //     for (i = 0; i < maxorder; ++i) {
    //         nbody_include[i] = i + 2;
    //     }
    // } else if (nbody_v.size() == maxorder) {
    //     for (i = 0; i < maxorder; ++i) {
    //         nbody_include[i] = boost::lexical_cast<int>(nbody_v[i]);
    //     }
    // }

    alm_core->interaction->maxorder = maxorder;
    alm_core->memory->allocate(alm_core->interaction->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i) {
        alm_core->interaction->nbody_include[i] = nbody_include[i];
    }
}

void ALM::set_cutoff_radii(const double * const * const * rcs)
{
    int i, j, k, nkd, maxorder;

    nkd = alm_core->system->nkd;
    maxorder = alm_core->interaction->maxorder;
    alm_core->memory->allocate(alm_core->interaction->rcs, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                alm_core->interaction->rcs[i][j][k] = rcs[i][j][k];
            }
        }
    }
}

ALMCore * ALM::get_alm_core()
{
    return alm_core;
}

int ALM::get_fc_length(const int fc_order)  // harmonic=2, ...
{
    int id, order, num_unique_elems, num_equiv_elems;
    Fcs *fcs;

    fcs = alm_core->fcs;
    id = 0;
    order = fc_order - 2;
    if (fcs->ndup[order].size() < 1) {return 0;}
    num_unique_elems = fcs->ndup[order].size();
    for (int iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
        num_equiv_elems = fcs->ndup[order][iuniq];
        id += num_equiv_elems;
    }
    return id;
}

void ALM::get_fc(double *fc_value,
                 int *elem_indices, // (len(fc_value), fc_order) is flatten.
                 const int fc_order) // harmonic=2, ...
{
    int j, k, ip, id;
    int order, num_unique_elems, num_equiv_elems, maxorder;
    double fc_elem, coef;
    Fcs *fcs;
    Fitting *fitting;

    fcs = alm_core->fcs;
    fitting = alm_core->fitting;
    maxorder = alm_core->interaction->maxorder;
    ip = 0;
    for (order = 0; order < maxorder; ++order) {
        if (fcs->ndup[order].size() < 1) {continue;} // 
        id = 0;
        num_unique_elems = fcs->ndup[order].size();
        for (int iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
            num_equiv_elems = fcs->ndup[order][iuniq];
            fc_elem = fitting->params[ip];
            for (j = 0; j < num_equiv_elems; ++j) {
                if (order == fc_order - 2) {
                    // coef is normally 1 or -1.
                    coef = fcs->fc_set[order][id].coef;
                    fc_value[id] = fc_elem * coef;
                    for (k = 0; k < fc_order; ++k) {
                        elem_indices[id * fc_order + k] = 
                            fcs->fc_set[order][id].elems[k];
                    }
                }
                ++id;
            }
            ++ip;
        }
        if (order >= fc_order - 2) {break;}
    }
}

void ALM::run_fitting()
{
    alm_core->constraint->setup();
    alm_core->fitting->fitmain();
}

void ALM::run_suggest()
{
    alm_core->displace->gen_displacement_pattern();
}


