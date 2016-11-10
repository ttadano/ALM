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
        InputParser(class ALM *, int, char **);
        ~InputParser();
        void parse_input(int, char **);

        std::string str_magmom;

    private:
        std::ifstream ifs_input;
        bool from_stdin;

        int locate_tag(std::string);
        void split_str_by_space(const std::string, std::vector<std::string> &);
        void parse_general_vars();
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
        void set_general_vars();
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
