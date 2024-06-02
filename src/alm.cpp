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
#include "optimize.h"
#include "cluster.h"
#include "patterndisp.h"
#include "symmetry.h"
#include "system.h"
#include "timer.h"
#include <iostream>
#include <string>
#include <Eigen/Sparse>


using namespace ALM_NS;

ALM::ALM()
{
    init_instances();
    verbosity = 1;
    structure_initialized = false;
    initialized_constraint_class = false;
    ofs_alm = nullptr;
    coutbuf = nullptr;
}

ALM::~ALM()
{
    delete files;
    delete system;
    delete cluster;
    delete fcs;
    delete symmetry;
    delete optimize;
    delete constraint;
    delete displace;
    delete timer;
    delete writer;
}

void ALM::init_instances()
{
    files = new Files();
    system = new System();
    cluster = new Cluster();
    fcs = new Fcs();
    symmetry = new Symmetry();
    optimize = new Optimize();
    constraint = new Constraint();
    displace = new Displace();
    timer = new Timer();
    writer = new Writer();
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
    files->set_prefix(prefix);
}

void ALM::set_print_symmetry(const int printsymmetry) const // PRINTSYM
{
    symmetry->set_print_symmetry(printsymmetry);
}

void ALM::set_datfile_train(const DispForceFile &dat_in) const
{
    files->set_datfile_train(dat_in);
}

void ALM::set_datfile_validation(const DispForceFile &dat_in) const
{
    files->set_datfile_validation(dat_in);
}

void ALM::set_symmetry_tolerance(const double tolerance) const // TOLERANCE
{
    symmetry->set_tolerance(tolerance);
}

void ALM::set_displacement_param(const bool trim_dispsign_for_evenfunc) const // TRIMEVEN
{
    displace->set_trim_dispsign_for_evenfunc(trim_dispsign_for_evenfunc);
}

void ALM::set_displacement_basis(const std::string str_disp_basis) const // DBASIS
{
    displace->set_disp_basis(str_disp_basis);
}

void ALM::set_periodicity(const int is_periodic[3]) const // PERIODIC
{
    system->set_periodicity(is_periodic);
}

void ALM::set_cell(const size_t nat,
                   const double lavec[3][3],
                   const double xcoord[][3],
                   const int kind[],
                   const std::string kdname[]) const
{
    system->set_supercell(lavec, nat, kind, xcoord);
    system->set_kdname(kdname);
}

void ALM::set_magnetic_params(const size_t nat,
                              const double (*magmom)[3], // MAGMOM
                              const bool lspin,
                              const int noncollinear, // NONCOLLINEAR
                              const int trev_sym_mag, // TREVSYM
                              const std::string str_magmom) const // MAGMOM
{
    system->set_spin_variables(nat,
                               lspin,
                               noncollinear,
                               trev_sym_mag,
                               magmom);
    system->set_str_magmom(str_magmom);
}

void ALM::set_u_train(const std::vector<std::vector<double>> &u) const
{
    optimize->set_u_train(u);
}

void ALM::set_f_train(const std::vector<std::vector<double>> &f) const
{
    optimize->set_f_train(f);
}

void ALM::set_validation_data(const std::vector<std::vector<double>> &u,
                              const std::vector<std::vector<double>> &f) const
{
    optimize->set_validation_data(u, f);
}

void ALM::set_optimizer_control(const OptimizerControl &optcontrol_in) const
{
    optimize->set_optimizer_control(optcontrol_in);
}

void ALM::set_constraint_mode(const int constraint_flag) const // ICONST
{
    constraint->set_constraint_mode(constraint_flag);
}

void ALM::set_algebraic_constraint(const int use_algebraic_flag) const // ICONST / 10
{
    constraint->set_constraint_algebraic(use_algebraic_flag);
}

void ALM::set_tolerance_constraint(const double tolerance_constraint) const // TOL_CONST
{
    constraint->set_tolerance_constraint(tolerance_constraint);
}

void ALM::set_rotation_axis(const std::string rotation_axis) const // ROTAXIS
{
    constraint->set_rotation_axis(rotation_axis);
}

void ALM::set_fc_file(const int order,
                      const std::string fc_file) const
{
    constraint->set_fc_file(order, fc_file);
}

void ALM::set_fc_fix(const int order,
                     const bool fc_fix) const
{
    if (order == 2) {
        constraint->set_fix_harmonic(fc_fix);
    }
    if (order == 3) {
        constraint->set_fix_cubic(fc_fix);
    }
}

bool ALM::ready_all_constraints() const
{
    return constraint->ready_all_constraints();
}

void ALM::set_forceconstants_to_fix(const std::vector<std::vector<int>> &intpair_fix,
                                    const std::vector<double> &values_fix) const
{
    constraint->set_forceconstants_to_fix(intpair_fix, values_fix);
}

void ALM::set_sparse_mode(const int sparse_mode) const // SPARSE
{
    auto optctrl = optimize->get_optimizer_control();
    optctrl.use_sparse_solver = sparse_mode;
    optimize->set_optimizer_control(optctrl);
}

void ALM::set_forceconstant_basis(const std::string preferred_basis) const // FCSYM_BASIS
{
    fcs->set_forceconstant_basis(preferred_basis);
}

std::string ALM::get_forceconstant_basis() const
{
    return fcs->get_forceconstant_basis();
}

void ALM::set_nmaxsave(const int nmaxsave) const // NMAXSAVE
{
    writer->set_output_maxorder(nmaxsave);
}

int ALM::get_nmaxsave() const
{
    return writer->get_output_maxorder();
}

void ALM::define(const int maxorder,
                 const size_t nkd,
                 const int *nbody_include,
                 const double *cutoff_radii) const
{
    // nkd = 0 means cutoff_radii undefined (hopefully nullptr).
    cluster->define(maxorder,
                    nkd,
                    nbody_include,
                    cutoff_radii);
}

OptimizerControl ALM::get_optimizer_control() const
{
    return optimize->get_optimizer_control();
}

std::vector<std::vector<double>> ALM::get_u_train() const
{
    return optimize->get_u_train();
}

std::vector<std::vector<double>> ALM::get_f_train() const
{
    return optimize->get_f_train();
}

size_t ALM::get_number_of_data() const
{
    return optimize->get_number_of_data();
}

size_t ALM::get_nrows_sensing_matrix() const
{
    return optimize->get_number_of_rows_sensing_matrix();
}

double ALM::get_cv_l1_alpha() const
{
    return optimize->get_cv_l1_alpha();
}

Cell ALM::get_supercell() const
{
    return system->get_supercell();
}

std::string *ALM::get_kdname() const
{
    return system->get_kdname();
}

Spin ALM::get_spin() const
{
    return system->get_spin();
}

std::string ALM::get_str_magmom() const
{
    return system->get_str_magmom();
}

double ***ALM::get_x_image() const
{
    return system->get_x_image();
}

int *ALM::get_periodicity() const
{
    return system->get_periodicity();
}

const std::vector<std::vector<int>> &ALM::get_atom_mapping_by_pure_translations() const
{
    return symmetry->get_map_p2s();
}

int ALM::get_maxorder() const
{
    return cluster->get_maxorder();
}

int *ALM::get_nbody_include() const
{
    return cluster->get_nbody_include();
}

size_t ALM::get_number_of_displacement_patterns(const int fc_order) const
// harmonic=1, ...
{
    const auto order = fc_order - 1;
    return displace->get_pattern_all(order).size();
}

void ALM::get_number_of_displaced_atoms(int *numbers,
                                        const int fc_order) const
// harmonic=1, ...
{
    const auto order = fc_order - 1;

    for (size_t i = 0; i < displace->get_pattern_all(order).size(); ++i) {
        numbers[i] = static_cast<int>(displace->get_pattern_all(order)[i].atoms.size());
    }
}

int ALM::get_displacement_patterns(int *atom_indices,
                                   double *disp_patterns,
                                   const int fc_order) const
// harmonic=1, ...
{
    const auto order = fc_order - 1;

    auto i_atom = 0;
    auto i_disp = 0;
    for (const auto &displacements: displace->get_pattern_all(order)) {
        for (size_t j = 0; j < displacements.atoms.size(); ++j) {
            atom_indices[i_atom] = displacements.atoms[j];
            ++i_atom;
            for (auto k = 0; k < 3; ++k) {
                disp_patterns[i_disp] = displacements.directions[3 * j + k];
                ++i_disp;
            }
        }
    }

    // 0:Cartesian or 1:Fractional. -1 means something wrong.
    if (displace->get_disp_basis()[0] == 'C') {
        return 0;
    }
    if (displace->get_disp_basis()[0] == 'F') {
        return 1;
    }
    return -1;
}

size_t ALM::get_number_of_fc_elements(const int fc_order) const
// harmonic=1, ...
{
    const auto order = fc_order - 1;

    if (fcs->get_nequiv()[order].empty()) { return 0; }
    size_t id = 0;
    const auto num_unique_elems = fcs->get_nequiv()[order].size();

    for (size_t iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
        const auto num_equiv_elems = fcs->get_nequiv()[order][iuniq];
        id += num_equiv_elems;
    }
    return id;
}

size_t ALM::get_number_of_irred_fc_elements(const int fc_order) // harmonic=1, ...
{
    // Returns the number of irreducible force constants for the given order.
    // The irreducible force constant means a set of independent force constants
    // reduced by using all available symmetry operations and
    // constraints for translational invariance. Rotational invariance is not considered.

    const auto order = fc_order - 1;
    if (!initialized_constraint_class) {
        constraint->setup(system,
                          fcs,
                          cluster,
                          symmetry,
                          get_optimizer_control().linear_model,
                          get_optimizer_control().mirror_image_conv,
                          verbosity,
                          timer);
        initialized_constraint_class = true;
    }
    if (!ready_all_constraints()) {
        constraint->update_constraint_matrix(system,
                                             symmetry,
                                             cluster,
                                             fcs,
                                             verbosity,
                                             get_optimizer_control().mirror_image_conv);
    }

    return constraint->get_index_bimap(order).size();
}

size_t ALM::get_number_of_fc_origin(const int fc_order,
                                    const int permutation) const
{
    if (fc_order <= 0) {
        std::cout << "fc_order must be larger than 0." << std::endl;
        exit(EXIT_FAILURE);
    }
    const auto maxorder = cluster->get_maxorder();
    if (fc_order > maxorder) {
        std::cout << "fc_order must not be larger than maxorder" << std::endl;
        exit(EXIT_FAILURE);
    }
    auto nfc_cart = fcs->get_nfc_cart(1);

    if (nfc_cart.size() < fc_order) {
        std::cout << "fc has not yet been computed or set." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (permutation) {
        return fcs->get_nfc_cart(1)[fc_order - 1];
    } else {
        return fcs->get_nfc_cart(0)[fc_order - 1];
    }
}

void ALM::get_fc_origin(double *fc_values,
                        int *elem_indices,  // (len(fc_values), fc_order + 1) is flatten.
                        const int fc_order, // harmonic=1, ...
                        const int permutation) const
{
    // Return a set of force constants Phi(i,j,k,...) where i is an atom
    // inside the primitive cell at origin.

    const auto maxorder = cluster->get_maxorder();
    if (fc_order > maxorder) {
        std::cout << "fc_order must not be larger than maxorder" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!fcs->get_fc_cart()) {
        std::cout << "fc has not yet been computed." << std::endl;
        exit(EXIT_FAILURE);
    }

    auto id = 0;

    if (permutation) {
        for (const auto &it: fcs->get_fc_cart()[fc_order - 1]) {
            fc_values[id] = it.fc_value;
            for (auto i = 0; i < fc_order + 1; ++i) {
                elem_indices[id * (fc_order + 1) + i] = it.flattenarray[i];
            }
            ++id;
        }
    } else {
        for (const auto &it: fcs->get_fc_cart()[fc_order - 1]) {
            if (it.is_ascending_order) {
                fc_values[id] = it.fc_value;
                for (auto i = 0; i < fc_order + 1; ++i) {
                    elem_indices[id * (fc_order + 1) + i] = it.flattenarray[i];
                }
                ++id;
            }
        }
    }
}


void ALM::get_fc_irreducible(double *fc_values,
                             int *elem_indices,  // (len(fc_values), fc_order + 1) is flatten.
                             const int fc_order) // harmonic=1, ...
{
    // Return an irreducible set of force constants.

    double fc_elem;

    const auto maxorder = cluster->get_maxorder();
    if (fc_order > maxorder) {
        std::cout << "fc_order must not be larger than maxorder" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!optimize->get_params()) {
        std::cout << "fc has not yet been computed." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!initialized_constraint_class) {
        constraint->setup(system,
                          fcs,
                          cluster,
                          symmetry,
                          get_optimizer_control().linear_model,
                          get_optimizer_control().mirror_image_conv,
                          verbosity,
                          timer);
        initialized_constraint_class = true;
    }
    if (!ready_all_constraints()) {
        constraint->update_constraint_matrix(system,
                                             symmetry,
                                             cluster,
                                             fcs,
                                             verbosity,
                                             get_optimizer_control().mirror_image_conv);
    }

    size_t ishift = 0;
    size_t inew, iold;

    for (auto order = 0; order < fc_order; ++order) {

        if (constraint->get_index_bimap(order).empty()) { continue; }

        if (order == fc_order - 1) {
            for (const auto &it: constraint->get_index_bimap(order)) {
                inew = it.left;
                iold = it.right + ishift;

                fc_elem = optimize->get_params()[iold];
                fc_values[inew] = fc_elem;
                for (auto i = 0; i < fc_order + 1; ++i) {
                    elem_indices[inew * (fc_order + 1) + i] =
                            fcs->get_fc_table()[order][it.right].elems[i];
                }
            }
        }
        ishift += fcs->get_nequiv()[order].size();
    }
}


void ALM::get_fc_all(double *fc_values,
                     int *elem_indices,  // (len(fc_values), fc_order + 1) is flatten.
                     const int fc_order, // harmonic=1, ...
                     const int permutation) const
{
    int i;
    const auto ntran = symmetry->get_ntran();
    const auto maxorder = cluster->get_maxorder();

    if (fc_order > maxorder) {
        std::cout << "fc_order must not be larger than maxorder" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!fcs->get_fc_cart()) {
        std::cout << "fc has not yet been computed." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<int> pair_tran(fc_order + 1);
    size_t id = 0;

    if (permutation) {
        for (const auto &it: fcs->get_fc_cart()[fc_order - 1]) {

            for (size_t itran = 0; itran < ntran; ++itran) {
                for (i = 0; i < fc_order + 1; ++i) {
                    pair_tran[i] = symmetry->get_map_sym()[it.atoms[i]][symmetry->get_symnum_tran()[itran]];
                }
                fc_values[id] = it.fc_value;
                for (i = 0; i < fc_order + 1; ++i) {
                    elem_indices[id * (fc_order + 1) + i] = 3 * pair_tran[i] + it.coords[i];
                }
                ++id;
            }
        }
    } else {
        for (const auto &it: fcs->get_fc_cart()[fc_order - 1]) {
            if (it.is_ascending_order) {
                for (size_t itran = 0; itran < ntran; ++itran) {
                    for (i = 0; i < fc_order + 1; ++i) {
                        pair_tran[i] = symmetry->get_map_sym()[it.atoms[i]][symmetry->get_symnum_tran()[itran]];
                    }
                    fc_values[id] = it.fc_value;
                    for (i = 0; i < fc_order + 1; ++i) {
                        elem_indices[id * (fc_order + 1) + i] = 3 * pair_tran[i] + it.coords[i];
                    }
                    ++id;
                }
            }
        }
    }
}

void ALM::get_fc_dependency_mat(const int fc_order,
                                int *elem_indices_irred,
                                int *elem_indices_origin,
                                double *matrix_out) const
{
    const auto maxorder = cluster->get_maxorder();

    if (fc_order > maxorder) {
        std::cout << "fc_order must not be larger than maxorder" << std::endl;
        exit(EXIT_FAILURE);
    }

    typedef Eigen::Triplet<double, size_t> T;
    std::vector<T> triplet_map, index_map_constraint, index_map_rotation;

    const auto order = fc_order - 1;
    const auto relation = constraint->get_index_bimap(order);
    const auto const_relate = constraint->get_const_relate(order);
    const auto nequiv = fcs->get_nequiv();
    const auto fcs_table = fcs->get_fc_table();

    if (relation.empty()) {
        std::cout << "All values are fixed, so the relation matrix cannot be generated" << std::endl;
        return;
    }

    for (const auto &it: relation) {
        size_t index_tmp = 0;
        for (auto i = 0; i < it.right; ++i) {
            index_tmp += nequiv[order][i];
        }
        for (auto i = 0; i < fc_order + 1; ++i) {
            elem_indices_irred[it.left * (fc_order + 1) + i]
                    = fcs_table[order][index_tmp].elems[i];
        }
        triplet_map.emplace_back(it.right, it.left, 1.0);
    }

    int irow_max = 0;
    int icol_max = 0;
    for (const auto &it: triplet_map) {
        if (it.row() > irow_max) irow_max = it.row();
        if (it.col() > icol_max) icol_max = it.col();
    }
    SpMat P_map(irow_max + 1, icol_max + 1);
    P_map.setFromTriplets(triplet_map.begin(), triplet_map.end());
    P_map.makeCompressed();

    std::set<size_t> set_first_index;
    // Equality constraint elements
    for (const auto &it: const_relate) {
        for (auto i = 0; i < it.alpha.size(); ++i) {
            index_map_constraint.emplace_back(it.p_index_target, it.p_index_orig[i], it.alpha[i]);
        }
        set_first_index.insert(it.p_index_target);
    }

    // Identity matrix elements
    const auto nrows_tmp = P_map.rows();
    for (auto i = 0; i < nrows_tmp; ++i) {
        if (set_first_index.find(i) == set_first_index.end()) {
            index_map_constraint.emplace_back(i, i, 1.0);
        }
    }

    irow_max = 0;
    icol_max = 0;
    for (const auto &it: index_map_constraint) {
        if (it.row() > irow_max) irow_max = it.row();
        if (it.col() > icol_max) icol_max = it.col();
    }
    SpMat M_constraint(irow_max + 1, icol_max + 1);
    M_constraint.setFromTriplets(index_map_constraint.begin(),
                                 index_map_constraint.end());
    M_constraint.makeCompressed();


    int ielems = 0;
    for (auto i = 0; i < nequiv[order].size(); ++i) {
        for (auto j = 0; j < nequiv[order][i]; ++j) {
            for (auto k = 0; k < fcs_table[order][ielems].elems.size(); ++k) {
                elem_indices_origin[ielems * (fc_order + 1) + k] = fcs_table[order][ielems].elems[k];
            }
            index_map_rotation.emplace_back(ielems, i, fcs_table[order][ielems].sign);
            ++ielems;
        }
    }

    SpMat M_rotation(ielems, nequiv[order].size());
    M_rotation.setFromTriplets(index_map_rotation.begin(), index_map_rotation.end());
    M_rotation.makeCompressed();

    SpMat M = M_rotation * M_constraint * P_map;
    Eigen::MatrixXd mat_dense = Eigen::MatrixXd(M);

    const auto ncol_irred = relation.size();

    //std::cout << "matrix size (C++) = " << mat_dense.rows() << "x" << mat_dense.cols() << '\n';

    for (int i = 0; i < mat_dense.rows(); ++i) {
//        std::cout << "row C++: ";
        for (int j = 0; j < ncol_irred; ++j) {
            matrix_out[i * ncol_irred + j] = mat_dense(i, j);
//            std::cout << std::setw(5) << mat_dense(i, j);
        }
//        std::cout << '\n';
    }

}

void ALM::set_fc(double *fc_in) const
{
    optimize->set_fcs_values(cluster->get_maxorder(),
                             fc_in,
                             fcs->get_nequiv(),
                             constraint);

    fcs->set_forceconstant_cartesian(cluster->get_maxorder(),
                                     optimize->get_params());
}

void ALM::set_fc_zero_threshold(const double threshold_in)
{
    fcs->set_fc_zero_threshold(threshold_in);
}

double ALM::get_fc_zero_threshold() const
{
    return fcs->get_fc_zero_threshold();
}

void ALM::get_matrix_elements(double *amat,
                              double *bvec)
{
    const auto maxorder = cluster->get_maxorder();
    double fnorm;

    std::vector<double> amat_vec;
    std::vector<double> bvec_vec;

    if (!initialized_constraint_class) {
        constraint->setup(system,
                          fcs,
                          cluster,
                          symmetry,
                          get_optimizer_control().linear_model,
                          get_optimizer_control().mirror_image_conv,
                          verbosity,
                          timer);
        initialized_constraint_class = true;
    }
    if (!ready_all_constraints()) {
        constraint->update_constraint_matrix(system,
                                             symmetry,
                                             cluster,
                                             fcs,
                                             verbosity,
                                             get_optimizer_control().mirror_image_conv);
    }

    optimize->get_matrix_elements_algebraic_constraint(maxorder,
                                                       amat_vec,
                                                       bvec_vec,
                                                       optimize->get_u_train(),
                                                       optimize->get_f_train(),
                                                       fnorm,
                                                       symmetry,
                                                       fcs,
                                                       constraint);
    // This may be inefficient.
    auto i = 0;
    for (const auto it: amat_vec) {
        amat[i++] = it;
    }
    i = 0;
    for (const auto it: bvec_vec) {
        bvec[i++] = it;
    }
    //amat = amat_vec.data();
    //bvec = bvec_vec.data();
}

int ALM::run_optimize()
{
    if (!structure_initialized) {
        std::cout << "initialize_structure must be called beforehand." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!initialized_constraint_class) {
        constraint->setup(system,
                          fcs,
                          cluster,
                          symmetry,
                          get_optimizer_control().linear_model,
                          get_optimizer_control().mirror_image_conv,
                          verbosity,
                          timer);
        initialized_constraint_class = true;
    }
    if (!ready_all_constraints()) {
        constraint->update_constraint_matrix(system,
                                             symmetry,
                                             cluster,
                                             fcs,
                                             verbosity,
                                             get_optimizer_control().mirror_image_conv);
    }

    const auto maxorder = cluster->get_maxorder();
    std::vector<std::string> str_order(maxorder);
    for (auto i = 0; i < maxorder; ++i) {
        str_order[i] = cluster->get_ordername(i);
    }
    const auto info = optimize->optimize_main(symmetry,
                                              constraint,
                                              fcs,
                                              maxorder,
                                              files->get_prefix(),
                                              str_order,
                                              verbosity,
                                              files->get_datfile_train(),
                                              files->get_datfile_validation(),
                                              writer->get_output_maxorder(),
                                              timer);
    return info;
}

void ALM::run_suggest()
{
    displace->gen_displacement_pattern(cluster,
                                       symmetry,
                                       fcs,
                                       constraint,
                                       system,
                                       verbosity);
}

void ALM::init_fc_table()
{
    // Initialization of structure information.
    // Perform initialization only once.

    if (!structure_initialized) {
        system->init(verbosity, timer);
        symmetry->init(system, verbosity, timer);
        structure_initialized = true;
    }

    // Build cluster & force constant table
    cluster->init(system,
                  symmetry,
                  get_optimizer_control().mirror_image_conv,
                  verbosity,
                  timer);
    fcs->init(cluster,
              symmetry,
              system->get_supercell(),
              verbosity,
              timer);

    // Switch off the initialized_constraint_class flag
    // because the force constants are updated
    // but corresponding constranits are not.
    initialized_constraint_class = false;
}

void ALM::save_fc(const std::string filename,
                  const std::string fcs_format,
                  const int maxorder_to_save) const
{
    writer->set_output_maxorder(maxorder_to_save);
    writer->set_filename_fcs(filename);
    writer->save_fcs_with_specific_format(fcs_format,
                                          system,
                                          symmetry,
                                          cluster,
                                          constraint,
                                          fcs,
                                          optimize,
                                          files,
                                          verbosity);
}

void ALM::set_fcs_save_flag(const std::string fcs_format, const int val) const
{
    writer->set_fcs_save_flag(fcs_format, val);
}

int ALM::get_fcs_save_flag(const std::string fcs_format) const
{
    return writer->get_fcs_save_flag(fcs_format);
}