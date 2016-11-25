/*
 alm.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <fstream>
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
#include "timer.h"
#include "patterndisp.h"

using namespace ALM_NS;

ALM::ALM()
{
    alm_core = new ALMCore();
    alm_core->create();
    verbose = true;
    ofs_alm = nullptr;
    coutbuf = nullptr;
}

ALM::~ALM()
{
    if (alm_core->system->kdname) {
        deallocate(alm_core->system->kdname);
    }
    alm_core->system->kdname = nullptr;
    if (alm_core->system->xcoord) {
        deallocate(alm_core->system->xcoord);
    }
    alm_core->system->xcoord = nullptr;
    if (alm_core->system->kd) {
        deallocate(alm_core->system->kd);
    }
    alm_core->system->kd = nullptr;
    if (alm_core->system->magmom) {
        deallocate(alm_core->system->magmom);
    }
    alm_core->system->magmom = nullptr;
    if (alm_core->interaction->nbody_include) {
        deallocate(alm_core->interaction->nbody_include);
    }
    alm_core->interaction->nbody_include = nullptr;
    if (alm_core->interaction->rcs) {
        deallocate(alm_core->interaction->rcs);
    }
    alm_core->interaction->rcs = nullptr;

    delete alm_core;
}

const void ALM::set_run_mode(const std::string mode)
{
    alm_core->mode = mode;
}

const void ALM::set_verbose(const bool verbose_in)
{
    verbose = verbose_in;
}

const void ALM::set_output_filename_prefix(const std::string prefix) // PREFIX
{
    alm_core->files->job_title = prefix;
}

const void ALM::set_is_print_symmetry(const int printsymmetry) // PRINTSYM
{
    alm_core->symmetry->printsymmetry = printsymmetry;
}

const void ALM::set_is_print_hessians(const bool print_hessian) // HESSIAN
{
    alm_core->files->print_hessian = print_hessian;
}

const void ALM::set_symmetry_param(const int nsym) // NSYM
{
    alm_core->symmetry->nsym = nsym;
}

const void ALM::set_symmetry_tolerance(const double tolerance) // TOLERANCE
{
    alm_core->symmetry->tolerance = tolerance;
}

const void ALM::set_displacement_param(const bool trim_dispsign_for_evenfunc) // TRIMEVEN
{
    alm_core->displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
}

const void ALM::set_displacement_basis(const std::string str_disp_basis) // DBASIS
{
    alm_core->displace->disp_basis = str_disp_basis;
}

const void ALM::set_periodicity(const int is_periodic[3]) // PERIODIC
{
    int i;
    for (i = 0; i < 3; ++i) {
        alm_core->interaction->is_periodic[i] = is_periodic[i];
    }
}

const void ALM::set_cell(const int nat,
                         const double lavec[3][3],
                         const double xcoord[][3],
                         const int kd[],
                         const std::string kdname[])
{
    int i, j, nkd;
    std::vector<int> nkd_vals(nat);
    bool kd_exist;

    nkd_vals[0] = kd[0];
    nkd = 1;
    for (i = 1; i < nat; ++i) {
        kd_exist = false;
        for (j = 0; j < nkd; ++j) {
            if (nkd_vals[j] == kd[i]) {
                kd_exist = true;
                break;
            }
        }
        if (!kd_exist) {
            nkd_vals[nkd] = kd[i];
            ++nkd;
        }
    }

    alm_core->system->nkd = nkd;
    alm_core->system->nat = nat;

    if (alm_core->system->xcoord) {
        deallocate(alm_core->system->xcoord);
    }
    allocate(alm_core->system->xcoord, nat, 3);

    if (alm_core->system->kd) {
        deallocate(alm_core->system->kd);
    }
    allocate(alm_core->system->kd, nat);

    if (alm_core->system->kdname) {
        deallocate(alm_core->system->kdname);
    }
    allocate(alm_core->system->kdname, nkd);

    if (alm_core->system->magmom) {
        deallocate(alm_core->system->magmom);
    }
    allocate(alm_core->system->magmom, nat, 3);

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
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            alm_core->system->magmom[i][j] = 0.0;
        }
    }
}

const void ALM::set_magnetic_params(const double *magmom, // MAGMOM
                                    const bool lspin, // MAGMOM
                                    const int noncollinear, // NONCOLLINEAR
                                    const int trev_sym_mag, // TREVSYM
                                    const std::string str_magmom) // MAGMOM
{
    int i, j, nat;

    nat = alm_core->system->nat;
    alm_core->system->lspin = lspin;
    alm_core->system->noncollinear = noncollinear;
    alm_core->system->str_magmom = str_magmom;
    alm_core->symmetry->trev_sym_mag = trev_sym_mag;

    if (alm_core->system->magmom) {
        deallocate(alm_core->system->magmom);
    }
    allocate(alm_core->system->magmom, nat, 3);

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            alm_core->system->magmom[i][j] = magmom[i * 3 + j];
        }
    }
}

const void ALM::set_displacement_and_force(const double *u_in,
                                           const double *f_in,
                                           const int nat,
                                           const int ndata_used)
{
    double **u;
    double **f;

    alm_core->system->ndata = ndata_used;
    alm_core->system->nstart = 1;
    alm_core->system->nend = ndata_used;

    allocate(u, ndata_used, 3 * nat);
    allocate(f, ndata_used, 3 * nat);

    for (int i = 0; i < ndata_used; i++) {
        for (int j = 0; j < 3 * nat; j++) {
            u[i][j] = u_in[i * nat * 3 + j];
            f[i][j] = f_in[i * nat * 3 + j];
        }
    }
    alm_core->fitting->set_displacement_and_force(u, f, nat, ndata_used);

    deallocate(u);
    deallocate(f);
}

const void ALM::set_fitting_constraint_type(const int constraint_flag) // ICONST
{
    alm_core->constraint->constraint_mode = constraint_flag;
}

const void ALM::set_fitting_constraint_rotation_axis
(const std::string rotation_axis) // ROTAXIS
{
    alm_core->constraint->rotation_axis = rotation_axis;
}

const void ALM::set_multiplier_option(const int multiply_data) // MULTDAT
{
    alm_core->symmetry->multiply_data = multiply_data;
}

const void ALM::set_fitting_filenames(const std::string dfile, // DFILE
                                      const std::string ffile) // FFILE
{
    alm_core->files->file_disp = dfile;
    alm_core->files->file_force = ffile;
}

const void ALM::set_norder(const int maxorder) // NORDER harmonic=1
{
    int i, j, k, nkd;

    alm_core->interaction->maxorder = maxorder;
    if (alm_core->interaction->nbody_include) {
        deallocate(alm_core->interaction->nbody_include);
    }
    allocate(alm_core->interaction->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i) {
        alm_core->interaction->nbody_include[i] = i + 2;
    }

    nkd = alm_core->system->nkd;
    if (alm_core->interaction->rcs) {
        deallocate(alm_core->interaction->rcs);
    }
    allocate(alm_core->interaction->rcs, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                alm_core->interaction->rcs[i][j][k] = -1.0;
            }
        }
    }
}

const void ALM::set_interaction_range(const int *nbody_include) // NBODY
{
    int maxorder = alm_core->interaction->maxorder;
    if (maxorder > 0) {
        if (alm_core->interaction->nbody_include) {
            deallocate(alm_core->interaction->nbody_include);
        }
        allocate(alm_core->interaction->nbody_include, maxorder);
        for (int i = 0; i < maxorder; ++i) {
            alm_core->interaction->nbody_include[i] = nbody_include[i];
        }
    }
}

const void ALM::set_cutoff_radii(const double *rcs)
{
    int i, j, k, nkd, maxorder, count;

    nkd = alm_core->system->nkd;
    maxorder = alm_core->interaction->maxorder;
    if (maxorder > 0) {
        if (alm_core->interaction->rcs) {
            deallocate(alm_core->interaction->rcs);
        }
        allocate(alm_core->interaction->rcs, maxorder, nkd, nkd);
    }

    count = 0;
    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                alm_core->interaction->rcs[i][j][k] = rcs[count];
                count++;
            }
        }
    }
}

ALMCore* ALM::get_alm_core()
{
    return alm_core;
}

const int ALM::get_atom_mapping_by_pure_translations(int *map_p2s)
{
    const int ntran = alm_core->symmetry->ntran;
    const int natmin = alm_core->symmetry->natmin;

    for (int i = 0; i < ntran; ++i) {
        for (int j = 0; j < natmin; ++j) {
            map_p2s[i * natmin + j] = alm_core->symmetry->map_p2s[j][i];
        }
    }
    return ntran;
}


const int ALM::get_number_of_displacement_patterns(const int fc_order) // harmonic=1, ...
{
    int order = fc_order - 1;

    return alm_core->displace->pattern_all[order].size();
}

const void ALM::get_numbers_of_displacements(int *numbers,
                                             const int fc_order) // harmonic=1, ...
{
    int order = fc_order - 1;

    for (int i = 0; i < alm_core->displace->pattern_all[order].size(); ++i) {
        numbers[i] = alm_core->displace->pattern_all[order][i].atoms.size();
    }
}

const int ALM::get_displacement_patterns(int *atom_indices,
                                         double *disp_patterns,
                                         const int fc_order) // harmonic=1, ...
{
    int i_atom, i_disp;
    AtomWithDirection *displacements;
    int order = fc_order - 1;

    i_atom = 0;
    i_disp = 0;
    for (int i = 0; i < alm_core->displace->pattern_all[order].size(); ++i) {
        displacements = &alm_core->displace->pattern_all[order][i];
        for (int j = 0; j < displacements->atoms.size(); ++j) {
            atom_indices[i_atom] = displacements->atoms[j];
            ++i_atom;
            for (int k = 0; k < 3; ++k) {
                disp_patterns[i_disp] = displacements->directions[3 * j + k];
                ++i_disp;
            }
        }
    }

    // 0:Cartesian or 1:Fractional. -1 means something wrong.
    if (alm_core->displace->disp_basis[0] == 'C') {
        return 0;
    } else if (alm_core->displace->disp_basis[0] == 'F') {
        return 1;
    } else {
        return -1;
    }
}

const int ALM::get_number_of_fc_elements(const int fc_order) // harmonic=1, ...
{
    int id, order, num_unique_elems, num_equiv_elems;
    Fcs *fcs;

    fcs = alm_core->fcs;
    id = 0;
    order = fc_order - 1;
    if (fcs->ndup[order].size() < 1) { return 0; }
    num_unique_elems = fcs->ndup[order].size();
    for (int iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
        num_equiv_elems = fcs->ndup[order][iuniq];
        id += num_equiv_elems;
    }
    return id;
}

const void ALM::get_fc(double *fc_values,
                       int *elem_indices, // (len(fc_values), fc_order + 1) is flatten.
                       const int fc_order) // harmonic=1, ...
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
    for (order = 0; order < fc_order; ++order) {
        if (fcs->ndup[order].size() < 1) { continue; }
        id = 0;
        num_unique_elems = fcs->ndup[order].size();
        for (int iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
            if (order == fc_order - 1) {
                fc_elem = fitting->params[ip];
                num_equiv_elems = fcs->ndup[order][iuniq];
                for (j = 0; j < num_equiv_elems; ++j) {
                    // coef is normally 1 or -1.
                    coef = fcs->fc_table[order][id].coef;
                    fc_values[id] = fc_elem * coef;
                    for (k = 0; k < fc_order + 1; ++k) {
                        elem_indices[id * (fc_order + 1) + k] =
                            fcs->fc_table[order][id].elems[k];
                    }
                    ++id;
                }
            }
            ++ip;
        }
    }
}

const void ALM::run()
{
    if (!verbose) {
        ofs_alm = new std::ofstream("alm.log", std::ofstream::out);
        coutbuf = std::cout.rdbuf();
        std::cout.rdbuf(ofs_alm->rdbuf());
    }

    alm_core->initialize();
    if (alm_core->mode == "fitting") {
        run_fitting();
    } else if (alm_core->mode == "suggest") {
        run_suggest();
    }

    if (!verbose) {
        ofs_alm->close();
        delete ofs_alm;
        ofs_alm = nullptr;
        std::cout.rdbuf(coutbuf);
        coutbuf = nullptr;
    }
}

const void ALM::run_fitting()
{
    alm_core->constraint->setup();
    alm_core->fitting->fitmain();
}

const void ALM::run_suggest()
{
    alm_core->displace->gen_displacement_pattern();
}


