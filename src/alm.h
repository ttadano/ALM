/*
 alm.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include "system.h"
#include "cluster.h"
#include "fcs.h"
#include "symmetry.h"
#include "optimize.h"
#include "constraint.h"
#include "files.h"
#include "patterndisp.h"
#include "timer.h"
#include "writer.h"

namespace ALM_NS {
class ALM {
public:
    ALM();

    ~ALM();

    class Cluster *cluster{};

    class Fcs *fcs{};

    class System *system{};

    class Symmetry *symmetry{};

    class Optimize *optimize{};

    class Constraint *constraint{};

    class Files *files{};

    class Displace *displace{};

    class Timer *timer{};

    class Writer *writer{};

    void set_verbosity(int verbosity_in);

    int get_verbosity() const;

    void set_output_filename_prefix(std::string prefix) const;

    void set_print_symmetry(int printsymmetry) const;

    void set_datfile_train(const DispForceFile &dat_in) const;

    void set_datfile_validation(const DispForceFile &dat_in) const;

    void set_symmetry_tolerance(double tolerance) const;

    void set_displacement_param(bool trim_dispsign_for_evenfunc) const;

    void set_displacement_basis(std::string str_disp_basis) const;

    void set_periodicity(const int is_periodic[3]) const;

    void set_cell(size_t nat,
                  const double lavec[3][3],
                  const double xcoord[][3],
                  const int kind[],
                  const std::string kdname[]) const;

    void set_magnetic_params(const size_t nat,
                             const double (*magmom)[3],
                             const bool lspin,
                             const int noncollinear,
                             const int trev_sym_mag,
                             const std::string str_magmom) const;

    void set_u_train(const std::vector<std::vector<double>> &u) const;

    void set_f_train(const std::vector<std::vector<double>> &f) const;

    void set_validation_data(const std::vector<std::vector<double>> &u,
                             const std::vector<std::vector<double>> &f) const;

    void set_optimizer_control(const OptimizerControl &optcontrol_in) const;

    void set_constraint_mode(const int constraint_flag) const;

    void set_algebraic_constraint(const int use_algebraic_flag) const;

    void set_tolerance_constraint(const double tolerance_constraint) const;

    void set_rotation_axis(const std::string rotation_axis) const;

    void set_fc_file(const int order, const std::string fc_file) const;

    void set_fc_fix(const int order, const bool fc_fix) const;

    bool ready_all_constraints() const;

    void set_forceconstants_to_fix(const std::vector<std::vector<int>> &intpair_fix,
                                   const std::vector<double> &values_fix) const;

    void set_sparse_mode(const int sparse_mode) const;

    void set_forceconstant_basis(const std::string preferred_basis) const;

    std::string get_forceconstant_basis() const;

    void set_nmaxsave(const int nmaxsave) const; // NMAXSAVE

    int get_nmaxsave() const;

    //void set_fitting_filenames(std::string dfile,
    //                           std::string ffile) const;
    void define(const int maxorder,
                const size_t nkd,
                const int *nbody_include,
                const double *cutoff_radii) const;

    //int get_ndata_used() const;
    OptimizerControl get_optimizer_control() const;

    std::vector<std::vector<double>> get_u_train() const;

    std::vector<std::vector<double>> get_f_train() const;

    size_t get_number_of_data() const;

    size_t get_nrows_sensing_matrix() const;

    double get_cv_l1_alpha() const;

    Cell get_supercell() const;

    std::string *get_kdname() const;

    Spin get_spin() const;

    void set_str_magmom(std::string);

    std::string get_str_magmom() const;

    double ***get_x_image() const;

    int *get_periodicity() const;

    const std::vector<std::vector<int>> &get_atom_mapping_by_pure_translations() const;

    int get_maxorder() const;

    int *get_nbody_include() const;

    size_t get_number_of_displacement_patterns(const int fc_order) const; // harmonic=1, ...
    void get_number_of_displaced_atoms(int *numbers,
                                       int fc_order) const; // harmonic=1, ...
    int get_displacement_patterns(int *atom_indices,
                                  double *disp_patterns,
                                  int fc_order) const;          // harmonic=1, ...
    size_t get_number_of_fc_elements(const int fc_order) const; // harmonic=1, ...
    size_t get_number_of_irred_fc_elements(const int fc_order); // harmonic=1, ...

    size_t get_number_of_fc_origin(const int fc_order, // harmonic = 1
                                   const int permutation) const;

    void get_fc_origin(double *fc_values,
                       int *elem_indices, // (len(fc_value), fc_order) is flatten.
                       int fc_order, // harmonic=1, ...
                       int permutation = 1) const;


    void get_fc_irreducible(double *fc_values,
                            int *elem_indices, // (len(fc_value), fc_order) is flatten.
                            int fc_order); // harmonic=1, ...


    void get_fc_all(double *fc_values,
                    int *elem_indices, // (len(fc_value), fc_order) is flatten.
                    int fc_order, // harmonic=1, ...
                    int permutation = 1) const;

    void get_fc_dependency_mat(const int fc_order,
                               int *index_elements_irred,
                               int *index_elements_origin,
                               double *matrix_out) const;

    void set_fc(double *fc_in) const;

    void set_fc_zero_threshold(const double threshold_in);

    double get_fc_zero_threshold() const;

    void get_matrix_elements(double *amat,
                             double *bvec);

    int run_optimize();

    void run_suggest();

    void init_fc_table();

    void save_fc(const std::string filename,
                 const std::string fc_format,
                 const int maxorder_to_save) const;

    void set_fcs_save_flag(const std::string fcs_format, const int val) const;

    int get_fcs_save_flag(const std::string fcs_format) const;

private:

    int verbosity;
    bool structure_initialized;
    bool initialized_constraint_class;
    std::ofstream *ofs_alm;
    std::streambuf *coutbuf;

    void init_instances();
};
}
