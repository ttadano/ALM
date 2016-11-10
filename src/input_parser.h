/*
 input_parser.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <fstream>
#include <string>
#include <map>
#include <vector>

namespace ALM_NS
{
    class InputParser: protected Pointers
    {
    public:
        InputParser(class ALMCore *);
        ~InputParser();
        void parse_input(int, char **);

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

        std::string str_magmom;

    private:
        std::ifstream ifs_input;
        bool from_stdin;

        int locate_tag(std::string);
        void split_str_by_space(const std::string, std::vector<std::string> &);
        void parse_general_vars();
        void parse_cell_parameter();
        void parse_interaction_vars();
        void parse_cutoff_radii();
        void parse_fitting_vars();
        void parse_atomic_positions();
        bool is_endof_entry(std::string);
        void get_var_dict(const std::string, std::map<std::string, std::string> &);

        template <typename T>
        void assign_val(T &, const std::string, std::map<std::string, std::string>);
    };
}
