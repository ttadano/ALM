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
//#include "alm.h"
#include "alm.h"
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
        void run(ALM *alm,
                 const int narg,
                 const char * const *arg);
        void parse_displacement_and_force(ALM *alm);
        void parse_displacement_and_force_files(Error *error,
                                                double **u,
                                                double **f,
                                                const int nat,
                                                const int ndata,
                                                const int nstart,
                                                const int nend,
                                                const std::string file_disp,
                                                const std::string file_force);
        std::string str_magmom;

    private:
        std::ifstream ifs_input;
        bool from_stdin;
        std::string *kdname;
        std::string mode;
        int maxorder;
        int nat;
        int nkd;

        void parse_input(ALM *alm);
        void parse_general_vars(ALM *alm);
        void parse_cell_parameter(ALM *alm);
        void parse_atomic_positions(ALM *alm);
        void parse_interaction_vars(ALM *alm);
        void parse_cutoff_radii(ALM *alm);
        void parse_fitting_vars(ALM *alm);
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
        void set_displacement_and_force(ALM *alm,
                                        const double * const *u,
                                        const double * const *f,
                                        const int nat,
                                        const int ndata_used,
                                        const int nmulti);
    };
}
