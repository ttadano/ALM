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
	void run(InputSetter *input,
		 const int narg,
		 const char * const *arg,
		 Error *error,
		 Memory *memory,
		 const std::string mode);
        std::string str_magmom;

    private:
        std::ifstream ifs_input;
        bool from_stdin;

	void parse_general_vars(InputSetter *input,
				Error *error,
				Memory *memory);
        void parse_cell_parameter(InputSetter *input,
				  Error *error);
        void parse_atomic_positions(InputSetter *input,
				    const int nat,
				    Error *error,
				    Memory *memory);
        int parse_interaction_vars(InputSetter *input,
				   Error *error,
				   Memory *memory);
        void parse_cutoff_radii(InputSetter *input,
				const int nkd,
				const int maxorder,
				const std::string *kdname,
				Error *error,
				Memory *memory);
        void parse_fitting_vars(InputSetter *input,
				Error *error);
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
    };
}
