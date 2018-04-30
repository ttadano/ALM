/*
 alm.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>

namespace ALM_NS
{
    class ALM
    {
    public:
        ALM();
        ~ALM();

        std::string mode;
        int verbosity;

        class InputSetter *input;
        class System *system;
        class Interaction *interaction;
        class Fcs *fcs;
        class Symmetry *symmetry;
        class Fitting *fitting;
        class Constraint *constraint;
        class Files *files;
        class Displace *displace;
        class Timer *timer;

        const void set_run_mode(std::string mode_in);
        const void set_verbose(bool verbose_in);
        const void set_output_filename_prefix(std::string prefix);
        const void set_is_print_symmetry(int is_printsymmetry);
        const void set_is_print_hessians(bool print_hessian);
        const void set_symmetry_param(int nsym);
        const void set_symmetry_tolerance(double tolerance);
        const void set_displacement_param(bool trim_dispsign_for_evenfunc);
        const void set_displacement_basis(std::string str_disp_basis);
        const void set_periodicity(const int is_periodic[3]);
        const void set_cell(int nat,
                            const double lavec[3][3],
                            const double xcoord[][3],
                            const int kd[],
                            const std::string kdname[]);
        const void set_magnetic_params(const double *magmom,
                                       bool lspin,
                                       int noncollinear,
                                       int trev_sym_mag,
                                       std::string str_magmom);
        const void set_displacement_and_force(const double *u_in,
                                              const double *f_in,
                                              int nat,
                                              int ndata_used);
        const int get_ndata_used();
        const void set_constraint_type(int constraint_flag);
        const void set_rotation_axis(std::string rotation_axis);
        const void set_sparse_mode(int sparse_mode);
        const void set_fitting_filenames(std::string dfile,
                                         std::string ffile);
        const void set_norder(int maxorder);
        const void set_nbody_include(const int *nbody_include);
        const void set_cutoff_radii(const double *rcs);

        const int get_atom_mapping_by_pure_translations(int *map_p2s);
        const int get_number_of_displacement_patterns(int fc_order); // harmonic=1, ...
        const void get_numbers_of_displacements(int *numbers,
                                                int fc_order); // harmonic=1, ...
        const int get_displacement_patterns(int *atom_indices,
                                            double *disp_patterns,
                                            int fc_order); // harmonic=1, ...
        const int get_number_of_fc_elements(int fc_order); // harmonic=2, ...
        const int get_number_of_irred_fc_elements(int fc_order); // harmonic=2, ...

        const void get_fc_origin(double *fc_value,
                                 int *elem_indices,
                                 // (len(fc_value), fc_order) is flatten.
                                 int fc_order); // harmonic=2, ...


        const void get_fc_irreducible(double *fc_value,
                                      int *elem_indices,
                                      // (len(fc_value), fc_order) is flatten.
                                      int fc_order); // harmonic=2, ...


        const void get_fc_all(double *fc_value,
                              int *elem_indices,
                              // (len(fc_value), fc_order) is flatten.
                              int fc_order); // harmonic=2, ...

        const void set_fc(double *fc_in);

        const void get_matrix_elements(int nat,
                                       int ndata_used,
                                       double *amat,
                                       double *bvec);
        const void generate_force_constant();
        const int optimize();
        const void run_suggest();
        const void run();

    private:

        bool structure_initialized;
        bool ready_to_fit;
        std::ofstream *ofs_alm;
        std::streambuf *coutbuf;
        void create();
        void initialize_structure(ALM *);
        void initialize_interaction(ALM *alm);
    };
}
