/*
 alm.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "alm.h"
#include "constraint.h"
#include "fcs.h"
#include "files.h"
#include "fitting.h"
#include "interaction.h"
#include "lasso.h"
#include "memory.h"
#include "patterndisp.h"
#include "symmetry.h"
#include "system.h"
#include "timer.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace ALM_NS;

ALM::ALM()
{
    create();
    verbosity = 1;
    structure_initialized = false;
    ready_to_fit = false;
    ofs_alm = nullptr;
    coutbuf = nullptr;
    mode = "suggest";
}

ALM::~ALM()
{
    if (interaction->nbody_include) {
        deallocate(interaction->nbody_include);
    }
    interaction->nbody_include = nullptr;
    if (interaction->cutoff_radii) {
        deallocate(interaction->cutoff_radii);
    }
    interaction->cutoff_radii = nullptr;

    delete files;
    delete system;
    delete interaction;
    delete fcs;
    delete symmetry;
    delete fitting;
    delete lasso;
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
    lasso = new Lasso();
    constraint = new Constraint();
    displace = new Displace();
    timer = new Timer();
}

void ALM::set_run_mode(const std::string mode_in)
{
    mode = mode_in;
}

std::string ALM::get_run_mode() const
{
    return mode;
}

void ALM::set_verbosity(const int verbosity_in)
{
    verbosity = verbosity_in;
}

int ALM::get_verbosity() const
{
    return verbosity;
}

void ALM::set_output_filename_prefix(const std::string prefix) const // PREFIX
{
    files->job_title = prefix;
}

void ALM::set_is_print_symmetry(const int printsymmetry) const // PRINTSYM
{
    symmetry->printsymmetry = printsymmetry;
}

void ALM::set_is_print_hessians(const bool print_hessian) const // HESSIAN
{
    files->print_hessian = print_hessian;
}

void ALM::set_symmetry_param(const int nsym) const // NSYM
{
    symmetry->nsym = nsym;
}

void ALM::set_symmetry_tolerance(const double tolerance) const // TOLERANCE
{
    symmetry->tolerance = tolerance;
}

void ALM::set_displacement_param(const bool trim_dispsign_for_evenfunc) const // TRIMEVEN
{
    displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
}

void ALM::set_displacement_basis(const std::string str_disp_basis) const // DBASIS
{
    displace->disp_basis = str_disp_basis;
}

void ALM::set_periodicity(const int is_periodic[3]) const // PERIODIC
{
    system->set_periodicity(is_periodic);
}

void ALM::set_cell(const int nat,
                   const double lavec[3][3],
                   const double xcoord[][3],
                   const int kd[],
                   const std::string kdname[]) const
{
    int i, j;
    std::vector<int> nkd_vals(nat);
    bool kd_exist;

    nkd_vals[0] = kd[0];
    int nkd = 1;
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

    // Generate the information of the supercell
    system->set_supercell(lavec, nat, nkd, kd, xcoord);
    system->set_kdname(kdname);
}

void ALM::set_magnetic_params(const double (*magmom)[3],
                              // MAGMOM
                              const bool lspin,
                              const int noncollinear,
                              // NONCOLLINEAR
                              const int trev_sym_mag,
                              // TREVSYM
                              const std::string str_magmom) const // MAGMOM
{
    system->set_spin_variables(lspin,
                               noncollinear,
                               trev_sym_mag,
                               magmom);
    system->set_str_magmom(str_magmom);
}

void ALM::set_displacement_and_force(const double *u_in,
                                     const double *f_in,
                                     const int nat,
                                     const int ndata_used) const
{
    double **u;
    double **f;

    fitting->ndata = ndata_used;
    fitting->nstart = 1;
    fitting->nend = ndata_used;

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

int ALM::get_ndata_used() const
{
    return fitting->get_ndata_used();
}

void ALM::set_constraint_type(const int constraint_flag) const // ICONST
{
    constraint->constraint_mode = constraint_flag;
}

void ALM::set_rotation_axis(const std::string rotation_axis) const // ROTAXIS
{
    constraint->rotation_axis = rotation_axis;
}

void ALM::set_sparse_mode(const int sparse_mode) const // SPARSE
{
    fitting->use_sparseQR = sparse_mode;
}

void ALM::set_fitting_filenames(const std::string dfile,
                                // DFILE
                                const std::string ffile) const // FFILE
{
    files->file_disp = dfile;
    files->file_force = ffile;
}

void ALM::set_norder(const int maxorder) const // NORDER harmonic=1
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

    nkd = system->get_supercell().number_of_elems;
    if (interaction->cutoff_radii) {
        deallocate(interaction->cutoff_radii);
    }
    allocate(interaction->cutoff_radii, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                interaction->cutoff_radii[i][j][k] = -1.0;
            }
        }
    }
}

void ALM::set_nbody_include(const int *nbody_include) const // NBODY
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

void ALM::set_cutoff_radii(const double *rcs) const
{
    const int nkd = system->get_supercell().number_of_elems;
    const auto maxorder = interaction->maxorder;

    if (maxorder > 0) {
        if (interaction->cutoff_radii) {
            deallocate(interaction->cutoff_radii);
        }
        allocate(interaction->cutoff_radii, maxorder, nkd, nkd);
    }

    auto count = 0;
    for (int i = 0; i < maxorder; ++i) {
        for (int j = 0; j < nkd; ++j) {
            for (int k = 0; k < nkd; ++k) {
                interaction->cutoff_radii[i][j][k] = rcs[count];
                count++;
            }
        }
    }
}


int ALM::get_atom_mapping_by_pure_translations(int *map_p2s) const
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


int ALM::get_number_of_displacement_patterns(const int fc_order) const
// harmonic=1, ...
{
    const auto order = fc_order - 1;
    return displace->pattern_all[order].size();
}

void ALM::get_number_of_displaced_atoms(int *numbers,
                                        const int fc_order) const
// harmonic=1, ...
{
    const auto order = fc_order - 1;

    for (int i = 0; i < displace->pattern_all[order].size(); ++i) {
        numbers[i] = displace->pattern_all[order][i].atoms.size();
    }
}

int ALM::get_displacement_patterns(int *atom_indices,
                                   double *disp_patterns,
                                   const int fc_order) const
// harmonic=1, ...
{
    AtomWithDirection *displacements;
    const auto order = fc_order - 1;

    auto i_atom = 0;
    auto i_disp = 0;
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
    }
    if (displace->disp_basis[0] == 'F') {
        return 1;
    }
    return -1;
}

int ALM::get_number_of_fc_elements(const int fc_order) const
// harmonic=1, ...
{
    const auto order = fc_order - 1;

    if (fcs->nequiv[order].empty()) { return 0; }
    auto id = 0;
    const int num_unique_elems = fcs->nequiv[order].size();

    for (int iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
        const auto num_equiv_elems = fcs->nequiv[order][iuniq];
        id += num_equiv_elems;
    }
    return id;
}

int ALM::get_number_of_irred_fc_elements(const int fc_order) // harmonic=1, ...
{
    // Returns the number of irreducible force constants for the given order.
    // The irreducible force constant means a set of independent force constants
    // reduced by using all available symmetry operations and
    // constraints for translational invariance. Rotational invariance is not considered.

    const auto order = fc_order - 1;
    if (!ready_to_fit) {
        constraint->setup(system,
                          fcs,
                          interaction,
                          symmetry,
                          mode,
                          verbosity,
                          timer);
        ready_to_fit = true;
    }
    return constraint->index_bimap[order].size();
}

void ALM::get_fc_origin(double *fc_values,
                        int *elem_indices,
                        // (len(fc_values), fc_order + 1) is flatten.
                        const int fc_order) const
// harmonic=1, ...
{
    // Return a set of force constants Phi(i,j,k,...) where i is an atom
    // inside the primitive cell at origin.

    int i;

    const auto maxorder = interaction->maxorder;
    if (fc_order > maxorder) {
        std::cout << "fc_order must not be larger than maxorder" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto ishift = 0;
    int ip;

    for (int order = 0; order < fc_order; ++order) {

        if (fcs->nequiv[order].empty()) { continue; }

        auto id = 0;

        if (order == fc_order - 1) {
            for (const auto &it : fcs->fc_table[order]) {

                ip = it.mother + ishift;
                fc_values[id] = fitting->params[ip] * it.sign;
                for (i = 0; i < fc_order + 1; ++i) {
                    elem_indices[id * (fc_order + 1) + i] = it.elems[i];
                }
                ++id;
            }
        }
        ishift += fcs->nequiv[order].size();
    }
}


void ALM::get_fc_irreducible(double *fc_values,
                             int *elem_indices,
                             // (len(fc_values), fc_order + 1) is flatten.
                             const int fc_order) // harmonic=1, ...
{
    // Return an irreducible set of force constants.

    int i;
    double fc_elem;

    const auto maxorder = interaction->maxorder;
    if (fc_order > maxorder) {
        std::cout << "fc_order must not be larger than maxorder" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!ready_to_fit) {
        constraint->setup(system,
                          fcs,
                          interaction,
                          symmetry,
                          mode,
                          verbosity,
                          timer);
        ready_to_fit = true;
    }

    auto ishift = 0;
    int inew, iold;

    for (int order = 0; order < fc_order; ++order) {

        if (constraint->index_bimap[order].empty()) { continue; }

        if (order == fc_order - 1) {
            for (const auto &it : constraint->index_bimap[order]) {
                inew = it.left;
                iold = it.right + ishift;

                fc_elem = fitting->params[iold];
                fc_values[inew] = fc_elem;
                for (i = 0; i < fc_order + 1; ++i) {
                    elem_indices[inew * (fc_order + 1) + i] =
                        fcs->fc_table[order][it.right].elems[i];
                }
            }
        }
        ishift += fcs->nequiv[order].size();
    }
}


void ALM::get_fc_all(double *fc_values,
                     int *elem_indices,
                     // (len(fc_values), fc_order + 1) is flatten.
                     const int fc_order) const
// harmonic=1, ...
{
    int i;
    double fc_elem;
    const int ntran = symmetry->ntran;

    const auto maxorder = interaction->maxorder;
    if (fc_order > maxorder) {
        std::cout << "fc_order must not be larger than maxorder" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto ishift = 0;
    int ip;
    std::vector<int> pair_tmp(fc_order + 1);
    std::vector<int> pair_tran(fc_order + 1);
    std::vector<int> xyz_tmp(fc_order + 1);

    for (int order = 0; order < fc_order; ++order) {

        if (fcs->nequiv[order].empty()) { continue; }

        auto id = 0;

        if (order == fc_order - 1) {
            for (const auto &it : fcs->fc_table[order]) {

                ip = it.mother + ishift;
                fc_elem = fitting->params[ip] * it.sign;

                for (i = 0; i < fc_order + 1; ++i) {
                    pair_tmp[i] = it.elems[i] / 3;
                    xyz_tmp[i] = it.elems[i] % 3;
                }

                for (int itran = 0; itran < ntran; ++itran) {
                    for (i = 0; i < fc_order + 1; ++i) {
                        pair_tran[i] = symmetry->map_sym[pair_tmp[i]][symmetry->symnum_tran[itran]];
                    }
                    fc_values[id] = fc_elem;
                    for (i = 0; i < fc_order + 1; ++i) {
                        elem_indices[id * (fc_order + 1) + i] =
                            3 * pair_tran[i] + xyz_tmp[i];
                    }
                    ++id;
                }
            }
        }

        ishift += fcs->nequiv[order].size();
    }
}

void ALM::set_fc(double *fc_in) const
{
    fitting->set_fcs_values(interaction->maxorder,
                            fc_in,
                            fcs->nequiv,
                            constraint);
}

void ALM::get_matrix_elements(const int ndata_used,
                              double *amat,
                              double *bvec) const
{
    const auto maxorder = interaction->maxorder;
    double fnorm;

    fitting->get_matrix_elements_algebraic_constraint(maxorder,
                                                      ndata_used,
                                                      amat,
                                                      bvec,
                                                      fnorm,
                                                      symmetry,
                                                      fcs,
                                                      constraint);
}


void ALM::generate_force_constant()
{
    initialize_structure();
    initialize_interaction();
}

void ALM::run()
{
    generate_force_constant();

    if (mode == "fitting") {
        optimize();
    } else if (mode == "suggest") {
        run_suggest();
    } else if (mode == "lasso") {
        optimize_lasso();
    }
}

int ALM::optimize()
{
    if (!structure_initialized) {
        std::cout << "initialize_structure must be called beforehand." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!ready_to_fit) {
        constraint->setup(system,
                          fcs,
                          interaction,
                          symmetry,
                          mode,
                          verbosity,
                          timer);
        ready_to_fit = true;
    }
    int info = fitting->fitmain(symmetry,
                                constraint,
                                fcs,
                                interaction->maxorder,
                                system->get_supercell().number_of_atoms,
                                verbosity,
                                files->file_disp,
                                files->file_force,
                                timer);
    return info;
}

void ALM::run_suggest() const
{
    displace->gen_displacement_pattern(interaction,
                                       symmetry,
                                       fcs,
                                       constraint,
                                       system,
                                       verbosity);
}

int ALM::optimize_lasso()
{
    if (!structure_initialized) {
        std::cout << "initialize_structure must be called beforehand." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!ready_to_fit) {
        constraint->setup(system,
                          fcs,
                          interaction,
                          symmetry,
                          mode,
                          verbosity,
                          timer);
        ready_to_fit = true;
    }
    lasso->lasso_main(symmetry,
                      interaction,
                      fcs,
                      constraint,
                      system->get_supercell().number_of_atoms,
                      files,
                      verbosity,
                      fitting,
                      timer);

    int info = 1;
    return info;
}


void ALM::initialize_structure()
{
    // Initialization of structure information.
    // Perform initialization only once.

    if (structure_initialized) return;
    system->init(verbosity, timer);
    files->init();
    symmetry->init(system, verbosity, timer);
    structure_initialized = true;
}

void ALM::initialize_interaction()
{
    // Build interaction & force constant table
    interaction->init(system,
                      symmetry,
                      verbosity,
                      timer);
    fcs->init(interaction,
              symmetry,
              system->get_supercell().number_of_atoms,
              verbosity,
              timer);

    // Switch off the ready flag because the force constants are updated
    // but corresponding constranits are not.
    ready_to_fit = false;
}
