/*
 input_parser.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm.h"
#include "input_setter.h"

#include <fstream>
#include <string>
#include <map>
#include <vector>


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

      //  void parse_displacement_and_force(ALM *alm) const;

        //void parse_displacement_and_force_files(double **u,
        //                                        double **f,
        //                                        const int nat_in,
        //                                        const int ndata,
        //                                        const int nstart,
        //                                        const int nend,
        //                                        const int skip_s,
        //                                        const int skip_e,
        //                                        const std::string file_disp,
        //                                        const std::string file_force) const;
        std::string str_magmom;

    private:
        std::ifstream ifs_input;
        bool from_stdin;
        std::string *kdname;
        std::string mode;
        int maxorder;
        int nat;
        int nkd;
        InputSetter *input_setter;

        void parse_input(ALM *alm);
        void parse_general_vars(ALM *alm);
        void parse_cell_parameter();
        void parse_atomic_positions(ALM *alm);
        void parse_interaction_vars();
        void parse_cutoff_radii();
        void parse_fitting_vars(ALM *alm);
        int locate_tag(const std::string);
        void split_str_by_space(const std::string,
                                std::vector<std::string> &) const;
        bool is_endof_entry(const std::string) const;
        void get_var_dict(const std::vector<std::string> &,
                          std::map<std::string,
                                   std::string> &);

        bool is_data_range_consistent(const DispForceFile &datfile_in) const;

        template <typename T>
        void assign_val(T &,
                        const std::string,
                        std::map<std::string, std::string>);

        void parse_displacement_and_force_files(std::vector<std::vector<double>> &u,
                                                std::vector<std::vector<double>> &f,
                                                DispForceFile &datfile_in) const;
    };
}
