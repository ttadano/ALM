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

        const void set_run_mode(const std::string mode);
        const void set_verbose(const bool verbose);
        const void set_output_filename_prefix(const std::string prefix);
        const void set_is_print_symmetry(const int is_printsymmetry);
        const void set_is_print_hessians(const bool print_hessian);
        const void set_symmetry_param(const int nsym);
        const void set_symmetry_tolerance(const double tolerance);
        const void set_displacement_param(const bool trim_dispsign_for_evenfunc);
        const void set_displacement_basis(const std::string str_disp_basis);
        const void set_periodicity(const int is_periodic[3]);
        const void set_cell(const int nat,
                            const double lavec[3][3],
                            const double xcoord[][3],
                            const int kd[],
                            const std::string kdname[]);
        const void set_magnetic_params(const double *magmom,
                                       const bool lspin,
                                       const int noncollinear,
                                       const int trev_sym_mag,
                                       const std::string str_magmom);
        const void set_displacement_and_force(const double *u_in,
                                              const double *f_in,
                                              const int nat,
                                              const int ndata_used);
        const void set_fitting_constraint_type(const int constraint_flag);
        const void set_fitting_constraint_rotation_axis
        (const std::string rotation_axis);
        const void set_multiplier_option(const int multiply_data);
        const void set_fitting_filenames(const std::string dfile,
                                         const std::string ffile);
        const void set_norder(const int maxorder);
        const void set_nbody_include(const int *nbody_include);
        const void set_cutoff_radii(const double *rcs);

        const int get_atom_mapping_by_pure_translations(int *map_p2s);
        const int get_number_of_displacement_patterns(const int fc_order); // harmonic=1, ...
        const void get_numbers_of_displacements(int *numbers,
                                                const int fc_order); // harmonic=1, ...
        const int get_displacement_patterns(int *atom_indices,
                                            double *disp_patterns,
                                            const int fc_order); // harmonic=1, ...
        const int get_number_of_fc_elements(const int fc_order); // harmonic=2, ...
        const int get_number_of_irred_fc_elements(const int fc_order); // harmonic=2, ...

        const void get_fc(double *fc_value,
                          int *elem_indices, // (len(fc_value), fc_order) is flatten.
                          const int fc_order); // harmonic=2, ...

        const void get_matrix_elements(const int nat, 
                                     const int ndata_used, 
                                     double *amat, 
                                     double *bvec);
        const void run();
        const void compute();

    private:

        bool verbose;
        std::ofstream *ofs_alm;
        std::streambuf *coutbuf;
        void create();
        void initialize(ALM *);
        void finalize();
        const void run_fitting(ALM *);
        const void run_suggest(ALM *);
    };
}
