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

        const void set_run_mode(const std::string mode);
        const void set_output_filename_prefix(const std::string prefix);
        const void set_is_print_symmetry(const int is_printsymmetry);
        const void set_is_print_hessians(const bool print_hessian);
        const void set_symmetry_params(const int nsym,
				       const double tolerance);
        const void set_displacement_params(const std::string str_disp_basis,
					   const bool trim_dispsign_for_evenfunc);
	const void set_displacement_basis(const std::string str_disp_basis);
        const void set_periodicity(const int is_periodic[3]);
        const void set_cell(const int nat,
			    const double lavec[3][3],
			    const double xcoord[][3],
			    const int kd[],
			    const std::string kdname[]);
        const void set_magnetic_params(const double* magmom,
				       const bool lspin,
				       const int noncollinear,
				       const int trev_sym_mag,
				       const std::string str_magmom);
        const void set_displacement_and_force(const double* u_in,
					      const double* f_in,
					      const int nat,
					      const int ndata_used);
        const void set_fitting_constraint(const int constraint_flag,
					  const std::string rotation_axis);
        const void set_multiplier_option(const int multiply_data);
        const void set_fitting_filenames(const std::string dfile,
					 const std::string ffile);
        const void set_norder(const int maxorder);
        const void set_interaction_range(const int *nbody_include);
	const void set_cutoff_radii(const double * rcs);
        ALMCore* get_alm_core();
        const int get_fc_length(const int fc_order); // harmonic=2, ...
        const void get_fc(double* fc_value,
			  int* elem_indices, // (len(fc_value), fc_order) is flatten.
			  const int fc_order); // harmonic=2, ...
        const void run();

    private:
        ALMCore* alm_core;
        const void run_fitting();
        const void run_suggest();
    };
}
