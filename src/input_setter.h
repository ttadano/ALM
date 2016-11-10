/*
 input_setter.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <string>

namespace ALM_NS
{
    class InputSetter: protected Pointers
    {
    public:
        InputSetter(class ALMCore *);
        ~InputSetter();
	void set_input(const std::string prefix,
		       const std::string mode,
		       const std::string str_disp_basis,
		       const std::string str_magmom,
		       const int nat,
		       const int nkd,
		       const int nsym,
		       const int is_printsymmetry,
		       const int is_periodic[3],
		       const bool trim_dispsign_for_evenfunc,
		       const bool lspin,
		       const bool print_hessian,
		       const int noncollinear,
		       const int trevsym,
		       const std::string *kdname,
		       const double * const *magmom,
		       const double tolerance,
		       const double a,
		       const double lavec_tmp[3][3],
		       const int maxorder,
		       const int *nbody_include,
		       const double * const * const * rcs,
		       const int ndata,
		       const int nstart,
		       const int nend,
		       const int nskip,
		       const int nboot,
		       const std::string dfile,
		       const std::string ffile,
		       const int multiply_data,
		       const int constraint_flag,
		       const std::string rotation_axis,
		       const std::string fc2_file,
		       const std::string fc3_file,
		       const bool fix_harmonic,
		       const bool fix_cubic,
		       const int *kd,
		       const double * const *xeq);
	void set_general_vars(const std::string prefix,
			      const std::string mode,
			      const std::string str_disp_basis,
			      const std::string str_magmom,
			      const int nat,
			      const int nkd,
			      const int nsym,
			      const int is_printsymmetry,
			      const int is_periodic[3],
			      const bool trim_dispsign_for_evenfunc,
			      const bool lspin,
			      const bool print_hessian,
			      const int noncollinear,
			      const int trevsym,
			      const std::string *kdname,
			      const double * const *magmom,
			      const double tolerance);
	void set_cell_parameter(const double a,
				const double lavec_tmp[3][3]);
	void set_interaction_vars(const int maxorder,
				  const int *nbody_include);
	void set_cutoff_radii(const int maxorder,
			      const int nkd,
			      const double * const * const * rcs);
	void set_fitting_vars(const int ndata,
			      const int nstart,
			      const int nend,
			      const int nskip,
			      const int nboot,
			      const std::string dfile,
			      const std::string ffile,
			      const int multiply_data,
			      const int constraint_flag,
			      const std::string rotation_axis,
			      const std::string fc2_file,
			      const std::string fc3_file,
			      const bool fix_harmonic,
			      const bool fix_cubic);
	void set_atomic_positions(const int nat,
				  const int *kd,
				  const double * const *xeq);
    };
}
