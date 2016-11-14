/*
 input_parser.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <fstream>
#include <string>
#include <map>
#include <vector>
#include "alm_core.h"
#include "input_setter.h"
#include "memory.h"
#include "error.h"

namespace ALM_NS
{
    class InputParser
    {
    public:
        InputParser();
        ~InputParser();
	void run(ALMCore *alm,
		 const int narg,
		 const char * const *arg);
	void parse_displacement_and_force(ALMCore *alm);
        std::string str_magmom;

    private:
        std::ifstream ifs_input;
        bool from_stdin;
	std::string *kdname;
	std::string mode;
	int maxorder;
	int nat;
	int nkd;

	void parse_input(ALMCore *alm);
	void parse_general_vars(ALMCore *alm);
	void parse_cell_parameter(ALMCore *alm);
        void parse_atomic_positions(ALMCore *alm);
        void parse_interaction_vars(ALMCore *alm);
        void parse_cutoff_radii(ALMCore *alm);
        void parse_fitting_vars(ALMCore *alm);
        int locate_tag(std::string);
        void split_str_by_space(const std::string, std::vector<std::string> &);
        bool is_endof_entry(std::string);
        void get_var_dict(const std::string,
			  std::map<std::string,
			  std::string> &,
			  Error *);

        template <typename T>
        void assign_val(T &,
			const std::string,
			std::map<std::string, std::string>,
			Error *);
	void data_multiplier(ALMCore *alm,
			     const int nat,
			     const int ndata,
			     const int nstart,
			     const int nend,
			     const int ndata_used,
			     const int multiply_data,
			     const std::string file_disp,
			     const std::string file_force);
    };
}
