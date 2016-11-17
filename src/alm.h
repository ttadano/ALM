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
    class ALMCore;

    class ALM
    {
    public:
        ALM();
        ~ALM();

        void set_run_mode(const std::string mode);
        void set_output_filename_prefix(const std::string prefix);
        void set_is_print_symmetry(const int is_printsymmetry);
        void set_is_print_hessians(const bool print_hessian);
        void set_symmetry_params(const int nsym,
                                 const double tolerance);
        void set_displacement_params(const std::string str_disp_basis,
                                     const bool trim_dispsign_for_evenfunc);
        void set_periodicity(const int is_periodic[3]);
        void set_cell(const int nat,
                      const int nkd,
                      const double lavec[3][3],
                      const double xcoord[][3],
                      const int kd[],
                      const std::string kdname[]);
        void set_magnetic_params(const double* const * magmom,
                                 const bool lspin,
                                 const int noncollinear,
                                 const int trev_sym_mag,
                                 const std::string str_magmom);
        void set_displacement_and_force(const double* u_in,
                                        const double* f_in,
                                        const int nat,
                                        const int ndata_used);
        void set_fitting_constraint(const int constraint_flag,
                                    const std::string rotation_axis);
        void set_multiplier_option(const int multiply_data);
        void set_fitting_filenames(const std::string dfile,
                                   const std::string ffile);
        void set_fitting_fc2_filename(const std::string fc2_file);
        void set_fitting_fc3_filename(const std::string fc3_file);
        void set_interaction_vars(const int maxorder,
                                  const int* nbody_include);
        void set_cutoff_radii(const double* const * const * rcs);

        ALMCore* get_alm_core();
        int get_fc_length(const int fc_order); // harmonic=2, ...
        void get_fc(double* fc_value,
                    int* elem_indices, // (len(fc_value), fc_order) is flatten.
                    const int fc_order); // harmonic=2, ...
        void run();

    private:
        ALMCore* alm_core;
        void run_fitting();
        void run_suggest();
    };
}
