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
//#include "alm.h"
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
    //alm = new ALMCore();
    //create();
    create();
    verbose = true;
    ofs_alm = nullptr;
    coutbuf = nullptr;
    mode = "suggest";
}

ALM::~ALM()
{
    if (system->kdname) {
        deallocate(system->kdname);
    }
    system->kdname = nullptr;
    if (system->xcoord) {
        deallocate(system->xcoord);
    }
    system->xcoord = nullptr;
    if (system->kd) {
        deallocate(system->kd);
    }
    system->kd = nullptr;
    if (system->magmom) {
        deallocate(system->magmom);
    }
    system->magmom = nullptr;
    if (interaction->nbody_include) {
        deallocate(interaction->nbody_include);
    }
    interaction->nbody_include = nullptr;
    if (interaction->rcs) {
        deallocate(interaction->rcs);
    }
    interaction->rcs = nullptr;

    delete files;
    delete system;
    delete interaction;
    delete fcs;
    delete symmetry;
    delete fitting;
    delete constraint;
    delete displace;
    delete timer;
}

void ALM::create()
{
    files = new Files();
    system = new System();
    interaction = new Interaction();
    fcs = new Fcs();
    symmetry = new Symmetry();
    fitting = new Fitting();
    constraint = new Constraint();
    displace = new Displace();
    timer = new Timer();
}

const void ALM::set_run_mode(const std::string mode_in)
{
    mode = mode_in;
}

const void ALM::set_verbose(const bool verbose_in)
{
    verbose = verbose_in;
}

const void ALM::set_output_filename_prefix(const std::string prefix) // PREFIX
{
    files->job_title = prefix;
}

const void ALM::set_is_print_symmetry(const int printsymmetry) // PRINTSYM
{
    symmetry->printsymmetry = printsymmetry;
}

const void ALM::set_is_print_hessians(const bool print_hessian) // HESSIAN
{
    files->print_hessian = print_hessian;
}

const void ALM::set_symmetry_param(const int nsym) // NSYM
{
    symmetry->nsym = nsym;
}

const void ALM::set_symmetry_tolerance(const double tolerance) // TOLERANCE
{
    symmetry->tolerance = tolerance;
}

const void ALM::set_displacement_param(const bool trim_dispsign_for_evenfunc) // TRIMEVEN
{
    displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
}

const void ALM::set_displacement_basis(const std::string str_disp_basis) // DBASIS
{
    displace->disp_basis = str_disp_basis;
}

const void ALM::set_periodicity(const int is_periodic[3]) // PERIODIC
{
    int i;
    for (i = 0; i < 3; ++i) {
        interaction->is_periodic[i] = is_periodic[i];
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

    system->nkd = nkd;
    system->nat = nat;

    if (system->xcoord) {
        deallocate(system->xcoord);
    }
    allocate(system->xcoord, nat, 3);

    if (system->kd) {
        deallocate(system->kd);
    }
    allocate(system->kd, nat);

    if (system->kdname) {
        deallocate(system->kdname);
    }
    allocate(system->kdname, nkd);

    if (system->magmom) {
        deallocate(system->magmom);
    }
    allocate(system->magmom, nat, 3);

    for (i = 0; i < nkd; ++i) {
        system->kdname[i] = kdname[i];
    }
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            system->lavec[i][j] = lavec[i][j];
        }
    }
    for (i = 0; i < nat; ++i) {
        system->kd[i] = kd[i];
        for (j = 0; j < 3; ++j) {
            system->xcoord[i][j] = xcoord[i][j];
        }
    }
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            system->magmom[i][j] = 0.0;
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

    nat = system->nat;
    system->lspin = lspin;
    system->noncollinear = noncollinear;
    system->str_magmom = str_magmom;
    symmetry->trev_sym_mag = trev_sym_mag;

    if (system->magmom) {
        deallocate(system->magmom);
    }
    allocate(system->magmom, nat, 3);

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            system->magmom[i][j] = magmom[i * 3 + j];
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

    system->ndata = ndata_used;
    system->nstart = 1;
    system->nend = ndata_used;

    allocate(u, ndata_used, 3 * nat);
    allocate(f, ndata_used, 3 * nat);

    for (int i = 0; i < ndata_used; i++) {
        for (int j = 0; j < 3 * nat; j++) {
            u[i][j] = u_in[i * nat * 3 + j];
            f[i][j] = f_in[i * nat * 3 + j];
        }
    }
    fitting->set_displacement_and_force(u, f, nat, ndata_used);

    deallocate(u);
    deallocate(f);
}

const void ALM::set_fitting_constraint_type(const int constraint_flag) // ICONST
{
    constraint->constraint_mode = constraint_flag;
}

const void ALM::set_fitting_constraint_rotation_axis
(const std::string rotation_axis) // ROTAXIS
{
    constraint->rotation_axis = rotation_axis;
}


const void ALM::set_fitting_filenames(const std::string dfile, // DFILE
                                      const std::string ffile) // FFILE
{
    files->file_disp = dfile;
    files->file_force = ffile;
}

const void ALM::set_norder(const int maxorder) // NORDER harmonic=1
{
    int i, j, k, nkd;

    interaction->maxorder = maxorder;
    if (interaction->nbody_include) {
        deallocate(interaction->nbody_include);
    }
    allocate(interaction->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i) {
        interaction->nbody_include[i] = i + 2;
    }

    nkd = system->nkd;
    if (interaction->rcs) {
        deallocate(interaction->rcs);
    }
    allocate(interaction->rcs, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                interaction->rcs[i][j][k] = -1.0;
            }
        }
    }
}

const void ALM::set_nbody_include(const int *nbody_include) // NBODY
{
    int maxorder = interaction->maxorder;
    if (maxorder > 0) {
        if (interaction->nbody_include) {
            deallocate(interaction->nbody_include);
        }
        allocate(interaction->nbody_include, maxorder);
        for (int i = 0; i < maxorder; ++i) {
            interaction->nbody_include[i] = nbody_include[i];
        }
    }
}

const void ALM::set_cutoff_radii(const double *rcs)
{
    int i, j, k, nkd, maxorder, count;

    nkd = system->nkd;
    maxorder = interaction->maxorder;
    if (maxorder > 0) {
        if (interaction->rcs) {
            deallocate(interaction->rcs);
        }
        allocate(interaction->rcs, maxorder, nkd, nkd);
    }

    count = 0;
    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                interaction->rcs[i][j][k] = rcs[count];
                count++;
            }
        }
    }
}

//ALMCore* ALM::get_alm()
//{
//    return alm;
//}

const int ALM::get_atom_mapping_by_pure_translations(int *map_p2s)
{
    const int ntran = symmetry->ntran;
    const int natmin = symmetry->nat_prim;

    for (int i = 0; i < ntran; ++i) {
        for (int j = 0; j < natmin; ++j) {
            map_p2s[i * natmin + j] = symmetry->map_p2s[j][i];
        }
    }
    return ntran;
}


const int ALM::get_number_of_displacement_patterns(const int fc_order) // harmonic=1, ...
{
    int order = fc_order - 1;

    return displace->pattern_all[order].size();
}

const void ALM::get_numbers_of_displacements(int *numbers,
                                             const int fc_order) // harmonic=1, ...
{
    int order = fc_order - 1;

    for (int i = 0; i < displace->pattern_all[order].size(); ++i) {
        numbers[i] = displace->pattern_all[order][i].atoms.size();
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
    for (int i = 0; i < displace->pattern_all[order].size(); ++i) {
        displacements = &displace->pattern_all[order][i];
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
    if (displace->disp_basis[0] == 'C') {
        return 0;
    } else if (displace->disp_basis[0] == 'F') {
        return 1;
    } else {
        return -1;
    }
}

const int ALM::get_number_of_fc_elements(const int fc_order) // harmonic=1, ...
{
    int id, order, num_unique_elems, num_equiv_elems;
    Fcs *fcs;

    fcs = fcs;
    id = 0;
    order = fc_order - 1;
    if (fcs->nequiv[order].size() < 1) { return 0; }
    num_unique_elems = fcs->nequiv[order].size();
    for (int iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
        num_equiv_elems = fcs->nequiv[order][iuniq];
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

    fcs = fcs;
    fitting = fitting;
    maxorder = interaction->maxorder;
    ip = 0;
    for (order = 0; order < fc_order; ++order) {
        if (fcs->nequiv[order].size() < 1) { continue; }
        id = 0;
        num_unique_elems = fcs->nequiv[order].size();
        for (int iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
            if (order == fc_order - 1) {
                fc_elem = fitting->params[ip];
                num_equiv_elems = fcs->nequiv[order][iuniq];
                for (j = 0; j < num_equiv_elems; ++j) {
                    // sign is normally 1 or -1.
                    coef = fcs->fc_table[order][id].sign;
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

const void ALM::run(ALM *alm)
{
    if (!verbose) {
        ofs_alm = new std::ofstream("alm.log", std::ofstream::out);
        coutbuf = std::cout.rdbuf();
        std::cout.rdbuf(ofs_alm->rdbuf());
    }

    initialize(alm);

    if (mode == "fitting") {
        run_fitting(alm);
    } else if (mode == "suggest") {
        run_suggest(alm);
    }

    if (!verbose) {
        ofs_alm->close();
        delete ofs_alm;
        ofs_alm = nullptr;
        std::cout.rdbuf(coutbuf);
        coutbuf = nullptr;
    }
}

const void ALM::run_fitting(ALM *alm)
{
    constraint->setup(alm);
    fitting->fitmain(alm);
}

const void ALM::run_suggest(ALM *alm)
{
    displace->gen_displacement_pattern(alm);
}


void ALM::initialize(ALM *alm)
{
    system->init(alm);
    files->init(alm);
    symmetry->init(alm);
    interaction->init(alm);
    fcs->init(alm);
}
