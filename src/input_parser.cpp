/*
 input_parser.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "input_parser.h"
#include "alm.h"
#include "error.h"
#include "files.h"
#include "optimize.h"
#include "input_setter.h"
#include "memory.h"
#include "system.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

InputParser::InputParser()
{
    input_setter = new InputSetter();
}

InputParser::~InputParser()
{
    delete input_setter;
}

void InputParser::run(ALM *alm,
                      const int narg,
                      const char * const *arg)
{
    if (narg == 1) {

        from_stdin = true;

    } else {

        from_stdin = false;

        ifs_input.open(arg[1], std::ios::in);
        if (!ifs_input) {
            std::cout << "No such file or directory: " << arg[1] << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    parse_input(alm);
}

void InputParser::parse_displacement_and_force(ALM *alm) const
{
    const int nat = alm->get_supercell().number_of_atoms;
    const auto ndata = alm->optimize->get_ndata();
    const auto nstart = alm->optimize->get_nstart();
    const auto nend = alm->optimize->get_nend();
    const auto skip_s = alm->optimize->get_skip_s();
    const auto skip_e = alm->optimize->get_skip_e();
    const auto ndata_used = nend - nstart + 1 - skip_e + skip_s;
    double **u;
    double **f;

    // Read displacement-force training data set from files
    const auto file_disp = alm->files->file_disp;
    const auto file_force = alm->files->file_force;

    allocate(u, ndata_used, 3 * nat);
    allocate(f, ndata_used, 3 * nat);
    parse_displacement_and_force_files(u, f, nat,
                                       ndata, nstart, nend,
                                       skip_s, skip_e,
                                       file_disp, file_force);
    alm->optimize->set_displacement_and_force(u, f, nat, ndata_used);
    deallocate(u);
    deallocate(f);
}

//
// This works independently from InputParser.
//
void InputParser::parse_displacement_and_force_files(double **u,
                                                     double **f,
                                                     const int nat_in,
                                                     const int ndata,
                                                     const int nstart,
                                                     const int nend,
                                                     const int skip_s,
                                                     const int skip_e,
                                                     const std::string file_disp,
                                                     const std::string file_force) const
{
    double u_in, f_in;
    std::ifstream ifs_disp, ifs_force;

    ifs_disp.open(file_disp.c_str(), std::ios::in);
    if (!ifs_disp) exit("openfiles", "cannot open disp file");
    ifs_force.open(file_force.c_str(), std::ios::in);
    if (!ifs_force) exit("openfiles", "cannot open force file");

    const unsigned int nreq = 3 * nat_in * ndata;

    std::vector<double> u_tmp(nreq), f_tmp(nreq);

    // Read displacements from DFILE

    unsigned int nline_u = 0;
    while (ifs_disp >> u_in) {
        u_tmp[nline_u++] = u_in;
        if (nline_u == nreq) break;
    }
    if (nline_u < nreq)
        exit("data_multiplier",
             "The number of lines in DFILE is too small for the given NDATA = ",
             ndata);

    // Read forces from FFILE

    unsigned int nline_f = 0;
    while (ifs_force >> f_in) {
        f_tmp[nline_f++] = f_in;
        if (nline_f == nreq) break;
    }
    if (nline_f < nreq)
        exit("data_multiplier",
             "The number of lines in FFILE is too small for the given NDATA = ",
             ndata);

    auto idata = 0;
    for (auto i = 0; i < ndata; ++i) {
        if (i < nstart - 1) continue;
        if (i >= skip_s && i < skip_e) continue;
        if (i > nend - 1) break;

        for (auto j = 0; j < nat_in; ++j) {
            for (auto k = 0; k < 3; ++k) {
                u[idata][3 * j + k] = u_tmp[3 * nat_in * i + 3 * j + k];
                f[idata][3 * j + k] = f_tmp[3 * nat_in * i + 3 * j + k];
            }
        }
        ++idata;
    }
    u_tmp.clear();
    f_tmp.clear();
    ifs_disp.close();
    ifs_force.close();
}

void InputParser::parse_input(ALM *alm)
{
    // The order of calling methods in this method is important.
    // Since following methods rely on variables those already
    // parsed.
    // Those below are set as the private class variables. See input_parser.h.
    //  std::string *kdname;
    //  std::string mode;
    //  int maxorder;
    //  int nat;
    //  int nkd;

    if (!locate_tag("&general")) {
        exit("parse_input",
             "&general entry not found in the input file");
    }

    // kdname is allocated in this method.
    parse_general_vars(alm);

    if (!locate_tag("&cell")) {
        exit("parse_input",
             "&cell entry not found in the input file");
    }
    parse_cell_parameter();

    if (!locate_tag("&position")) {
        exit("parse_input",
             "&position entry not found in the input file");
    }
    parse_atomic_positions(alm);
    input_setter->set_geometric_structure(alm);

    if (!locate_tag("&interaction")) {
        exit("parse_input",
             "&interaction entry not found in the input file");
    }
    parse_interaction_vars();

    if (!locate_tag("&cutoff")) {
        exit("parse_input",
             "&cutoff entry not found in the input file");
    }
    parse_cutoff_radii();
    input_setter->define(alm);

    if (mode == "optimize") {
        if (!locate_tag("&optimize")) {
            if (!locate_tag("&fitting")) {
                exit("parse_input",
                    "&optimize entry not found in the input file");
            } else {
                warn("parse_input", "&fitting field is deprecated. Please use &optimize instead.");
            }
        }
        parse_fitting_vars(alm);
    }

    deallocate(kdname);
}

void InputParser::parse_general_vars(ALM *alm)
{
    int i, j;
    std::string str_tmp, str_disp_basis;
    int printsymmetry, is_periodic[3];
    int icount, ncount;
    auto trim_dispsign_for_evenfunc = true;
    bool print_hessian;
    int noncollinear, trevsym;
    double **magmom, magmag;
    double tolerance;
    double tolerance_constraint;
    int verbosity;

    std::vector<std::string> kdname_v, periodic_v, magmom_v, str_split;
    std::string str_allowed_list =
        "PREFIX MODE NAT NKD KD PERIODIC PRINTSYM TOLERANCE DBASIS TRIMEVEN\
                                   MAGMOM NONCOLLINEAR TREVSYM HESSIAN TOL_CONST VERBOSITY";
    std::string str_no_defaults = "PREFIX MODE NAT NKD KD";
    std::vector<std::string> no_defaults;
    std::map<std::string, std::string> general_var_dict;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(str_allowed_list, general_var_dict);

    boost::split(no_defaults, str_no_defaults, boost::is_space());

    for (const auto &it : no_defaults) {
        if (general_var_dict.find(it) == general_var_dict.end()) {
            exit("parse_general_vars",
                 "The following variable is not found in &general input region: ",
                 it.c_str());
        }
    }

    const auto prefix = general_var_dict["PREFIX"];
    mode = general_var_dict["MODE"];

    std::transform(mode.begin(), mode.end(), mode.begin(), tolower);
    if (mode != "fitting" && mode != "suggest" && mode != "lasso" && mode != "optimize") {
        exit("parse_general_vars", "Invalid MODE variable");
    }
    if (mode == "fitting") {
        mode = "optimize";
        warn("parse_general_vars", "MODE = fitting is deprecated. Please use MODE = optimize instead.");
    }
    if (mode == "lasso") {
        mode = "optimize";
        warn("parse_general_vars", 
            "MODE = lasso is deprecated. Please use MODE = optimize instead with OPTIMIZER = enet option in the &optimize field.");
    }

    assign_val(nat, "NAT", general_var_dict);
    assign_val(nkd, "NKD", general_var_dict);

    if (general_var_dict["VERBOSITY"].empty()) {
        verbosity = 1;
    } else {
        assign_val(verbosity, "VERBOSITY", general_var_dict);
    }

    if (general_var_dict["PRINTSYM"].empty()) {
        printsymmetry = 0;
    } else {
        assign_val(printsymmetry, "PRINTSYM", general_var_dict);
    }

    split_str_by_space(general_var_dict["KD"], kdname_v);

    if (kdname_v.size() != nkd) {
        exit("parse_general_vars",
             "The number of entries for KD is inconsistent with NKD");
    } else {
        allocate(kdname, nkd);
        for (i = 0; i < nkd; ++i) {
            kdname[i] = kdname_v[i];
        }
    }

    split_str_by_space(general_var_dict["PERIODIC"], periodic_v);

    if (periodic_v.empty()) {
        for (i = 0; i < 3; ++i) {
            is_periodic[i] = 1;
        }
    } else if (periodic_v.size() == 3) {
        for (i = 0; i < 3; ++i) {
            try {
                is_periodic[i] = boost::lexical_cast<int>(periodic_v[i]);
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                exit("parse_general_vars",
                     "The PERIODIC tag must be a set of integers.");
            }
        }
    } else {
        exit("parse_general_vars",
             "Invalid number of entries for PERIODIC");
    }

    if (general_var_dict["TOLERANCE"].empty()) {
        tolerance = 1.0e-6;
    } else {
        assign_val(tolerance, "TOLERANCE", general_var_dict);
    }
    if (general_var_dict["TOL_CONST"].empty()) {
        tolerance_constraint = 1.0e-6;
    } else {
        assign_val(tolerance_constraint, "TOL_CONST", general_var_dict);
    }

    // Convert MAGMOM input to array
    allocate(magmom, nat, 3);
    auto lspin = false;

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            magmom[i][j] = 0.0;
        }
    }

    if (general_var_dict["NONCOLLINEAR"].empty()) {
        noncollinear = 0;
    } else {
        assign_val(noncollinear, "NONCOLLINEAR", general_var_dict);
    }
    if (general_var_dict["TREVSYM"].empty()) {
        trevsym = 1;
    } else {
        assign_val(trevsym, "TREVSYM", general_var_dict);
    }
    if (general_var_dict["HESSIAN"].empty()) {
        print_hessian = false;
    } else {
        assign_val(print_hessian, "HESSIAN", general_var_dict);
    }

    if (!general_var_dict["MAGMOM"].empty()) {
        lspin = true;

        if (noncollinear) {
            icount = 0;
            split_str_by_space(general_var_dict["MAGMOM"], magmom_v);
            for (std::vector<std::string>::const_iterator it = magmom_v.begin();
                 it != magmom_v.end(); ++it) {
                if ((*it).find('*') != std::string::npos) {
                    exit("parse_general_vars",
                         "Wild card '*' is not supported when NONCOLLINEAR = 1.");
                } else {
                    magmag = boost::lexical_cast<double>((*it));
                    if (icount / 3 >= nat) {
                        exit("parse_general_vars", "Too many entries for MAGMOM.");
                    }
                    magmom[icount / 3][icount % 3] = magmag;
                    ++icount;
                }
            }

            if (icount != 3 * nat) {
                exit("parse_general_vars",
                     "Number of entries for MAGMOM must be 3*NAT when NONCOLLINEAR = 1.");
            }
        } else {
            icount = 0;
            split_str_by_space(general_var_dict["MAGMOM"], magmom_v);
            for (std::vector<std::string>::const_iterator it = magmom_v.begin();
                 it != magmom_v.end(); ++it) {

                if ((*it).find('*') != std::string::npos) {
                    if ((*it) == "*") {
                        exit("parse_general_vars",
                             "Please place '*' without space for the MAGMOM-tag.");
                    }
                    boost::split(str_split, (*it), boost::is_any_of("*"));
                    if (str_split.size() != 2) {
                        exit("parse_general_vars",
                             "Invalid format for the MAGMOM-tag.");
                    } else {
                        if (str_split[0].empty() || str_split[1].empty()) {
                            exit("parse_general_vars",
                                 "Please place '*' without space for the MAGMOM-tag.");
                        }
                        try {
                            magmag = boost::lexical_cast<double>(str_split[1]);
                            ncount = static_cast<int>(boost::lexical_cast<double>(str_split[0]));
                        }
                        catch (std::exception &e) {
                            exit("parse_general_vars", "Bad format for MAGMOM.");
                        }

                        for (i = icount; i < icount + ncount; ++i) {
                            magmom[i][2] = magmag;
                        }
                        icount += ncount;
                    }

                } else {
                    magmag = boost::lexical_cast<double>((*it));
                    if (icount == nat) {
                        icount = 0;
                        break;
                    }
                    magmom[icount++][2] = magmag;
                }
            }
            if (icount != nat) {
                exit("parse_general_vars",
                     "Number of entries for MAGMOM must be NAT.");
            }
        }
    }

    if (mode == "suggest") {
        if (general_var_dict["DBASIS"].empty()) {
            str_disp_basis = "Cart";
        } else {
            str_disp_basis = general_var_dict["DBASIS"];
        }
        std::transform(str_disp_basis.begin(), str_disp_basis.end(),
                       str_disp_basis.begin(), toupper);
        if ((str_disp_basis[0] != 'C') && (str_disp_basis[0] != 'F')) {
            exit("parse_general_vars", "Invalid DBASIS");
        }

        if (!general_var_dict["TRIMEVEN"].empty()) {
            assign_val(trim_dispsign_for_evenfunc,
                       "TRIMEVEN", general_var_dict);
        }

    }

    input_setter->set_general_vars(alm,
                                   prefix,
                                   mode,
                                   verbosity,
                                   str_disp_basis,
                                   general_var_dict["MAGMOM"],
                                   nat,
                                   nkd,
                                   printsymmetry,
                                   is_periodic,
                                   trim_dispsign_for_evenfunc,
                                   lspin,
                                   print_hessian,
                                   noncollinear,
                                   trevsym,
                                   kdname,
                                   magmom,
                                   tolerance,
                                   tolerance_constraint);
    allocate(magmom, nat, 3);

    kdname_v.clear();
    periodic_v.clear();
    no_defaults.clear();
    general_var_dict.clear();
}

void InputParser::parse_cell_parameter()
{
    auto a = 0.0;
    double lavec_tmp[3][3];
    std::string line;
    std::string line_wo_comment, line_tmp;
    std::vector<std::string> line_vec, line_split;
    std::string::size_type pos_first_comment_tag;

    line_vec.clear();

    if (from_stdin) {

        while (std::getline(std::cin, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\n\r "));

            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            line_vec.push_back(line_wo_comment);
        }

    } else {
        while (std::getline(ifs_input, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\n\r "));

            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            line_vec.push_back(line_wo_comment);
        }
    }

    if (line_vec.size() != 4) {
        exit("parse_cell_parameter",
             "Too few or too much lines for the &cell field.\n \
                                            The number of valid lines for the &cell field should be 4.");
    }

    for (auto i = 0; i < 4; ++i) {

        line = line_vec[i];
        boost::split(line_split, line, boost::is_any_of("\t "), boost::token_compress_on);

        if (i == 0) {
            // Lattice factor a
            if (line_split.size() == 1) {
                a = boost::lexical_cast<double>(line_split[0]);
            } else {
                exit("parse_cell_parameter",
                     "Unacceptable format for &cell field.");
            }

        } else {
            // Lattice vectors a1, a2, a3
            if (line_split.size() == 3) {
                for (auto j = 0; j < 3; ++j) {
                    lavec_tmp[j][i - 1] = boost::lexical_cast<double>(line_split[j]);
                }
            } else {
                exit("parse_cell_parameter",
                     "Unacceptable format for &cell field.");
            }
        }
    }

    input_setter->set_cell_parameter(a, lavec_tmp);
}


void InputParser::parse_interaction_vars()
{
    int i;
    int *nbody_include;

    std::vector<std::string> nbody_v;
    std::string str_allowed_list = "NORDER NBODY";
    std::string str_no_defaults = "NORDER";
    std::vector<std::string> no_defaults;
    std::map<std::string, std::string> interaction_var_dict;


    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(str_allowed_list, interaction_var_dict);

    boost::split(no_defaults, str_no_defaults, boost::is_space());

    for (const auto &it : no_defaults) {
        if (interaction_var_dict.find(it) == interaction_var_dict.end()) {
            exit("parse_interaction_vars",
                 "The following variable is not found in &interaction input region: ",
                 it.c_str());
        }
    }

    assign_val(maxorder, "NORDER", interaction_var_dict);
    if (maxorder < 1)
        exit("parse_interaction_vars",
             "maxorder has to be a positive integer");

    allocate(nbody_include, maxorder);

    boost::split(nbody_v, interaction_var_dict["NBODY"], boost::is_space());

    if (nbody_v[0].empty()) {
        for (i = 0; i < maxorder; ++i) {
            nbody_include[i] = i + 2;
        }
    } else if (nbody_v.size() == maxorder) {
        for (i = 0; i < maxorder; ++i) {
            try {
                nbody_include[i] = boost::lexical_cast<int>(nbody_v[i]);
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                exit("parse_interaction_vars",
                     "NBODY must be an integer.");
            }
        }
    } else {
        exit("parse_interaction_vars",
             "The number of entry of NBODY has to be equal to NORDER");
    }

    if (nbody_include[0] != 2) {
        warn("parce_input",
             "Harmonic interaction is always 2 body (except on-site 1 body)");
    }


    input_setter->set_interaction_vars(maxorder, nbody_include);

    deallocate(nbody_include);

    nbody_v.clear();
    no_defaults.clear();
}


void InputParser::parse_fitting_vars(ALM *alm)
{
    int ndata, nstart, nend;
    int constraint_flag;
    auto flag_sparse = 0;
    std::string rotation_axis;
    std::string str_allowed_list =
        "NDATA NSTART NEND DFILE FFILE ICONST ROTAXIS FC2XML FC3XML SPARSE \
                                   LASSO_DNORM LASSO_ALPHA LASSO_MAXITER LASSO_TOL LASSO_CV LASSO_CVSET \
                                   LASSO_FREQ LASSO_MAXALPHA LASSO_MINALPHA LASSO_NALPHA \
                                   NDATA_TEST DFILE_TEST FFILE_TEST NSTART_TEST NEND_TEST SKIP STANDARDIZE \
                                   SOLUTION_PATH DEBIAS_OLS L1_RATIO OPTIMIZER";
    std::string str_no_defaults = "NDATA DFILE FFILE";
    std::vector<std::string> no_defaults;
    std::map<std::string, std::string> fitting_var_dict;

    int ndata_test, nstart_test, nend_test;
    int skip_s = 0;
    int skip_e = 0;
    std::string str_skip, dfile_test, ffile_test;
    OptimizerControl optcontrol;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(str_allowed_list, fitting_var_dict);

    boost::split(no_defaults, str_no_defaults, boost::is_space());

    for (const auto &it : no_defaults) {
        if (fitting_var_dict.find(it) == fitting_var_dict.end()) {
            exit("parse_fitting_vars",
                 "The following variable is not found in &fitting input region: ",
                 it.c_str());
        }
    }

    if (!fitting_var_dict["OPTIMIZER"].empty()) {
        auto str_optimizer = fitting_var_dict["OPTIMIZER"];
        boost::to_lower(str_optimizer);

        if (str_optimizer == "ols" || str_optimizer == "ls" || str_optimizer == "least-squares") {
            optcontrol.optimizer = 1;
        } else if (str_optimizer == "enet" || str_optimizer == "elastic-net") {
            optcontrol.optimizer = 2;
        } else {
            exit("parse_fitting_vars", "Invalid OPTIMIZER-tag");
        }
    }

    assign_val(ndata, "NDATA", fitting_var_dict);

    if (fitting_var_dict["NSTART"].empty()) {
        nstart = 1;
    } else {
        assign_val(nstart, "NSTART", fitting_var_dict);
    }
    if (fitting_var_dict["NEND"].empty()) {
        nend = ndata;
    } else {
        assign_val(nend, "NEND", fitting_var_dict);
    }

    if (!fitting_var_dict["SKIP"].empty()) {
        std::vector<std::string> str_entry;
        assign_val(str_skip, "SKIP", fitting_var_dict);
        boost::split(str_entry, str_skip, boost::is_any_of("-"));
        if (str_entry.size() == 1) {
            skip_s = boost::lexical_cast<int>(str_entry[0]) - 1;
            skip_e = skip_s + 1;
        } else if (str_entry.size() == 2) {
            skip_s = boost::lexical_cast<int>(str_entry[0]) - 1;
            skip_e = boost::lexical_cast<int>(str_entry[1]);
        } else {
            exit("parse_fitting_vars", "Invalid format for the SKIP-tag.");
        }
    }


    if (ndata <= 0 || nstart <= 0 || nend <= 0
        || nstart > ndata || nend > ndata || nstart > nend) {
        exit("parce_fitting_vars",
             "ndata, nstart, nend are not consistent with each other");
    }

    auto dfile = fitting_var_dict["DFILE"];
    auto ffile = fitting_var_dict["FFILE"];

    if (fitting_var_dict["ICONST"].empty()) {
        constraint_flag = 1;
    } else {
        assign_val(constraint_flag, "ICONST", fitting_var_dict);
    }

    auto fc2_file = fitting_var_dict["FC2XML"];
    auto fc3_file = fitting_var_dict["FC3XML"];
    bool fix_harmonic = !fc2_file.empty();
    bool fix_cubic = !fc3_file.empty();

    if (constraint_flag % 10 >= 2) {
        rotation_axis = fitting_var_dict["ROTAXIS"];
        if (rotation_axis.empty()) {
            exit("parse_fitting_vars",
                 "ROTAXIS has to be given when ICONST=2 or 3");
        }
    }

    if (!fitting_var_dict["SPARSE"].empty()) {
        assign_val(flag_sparse, "SPARSE", fitting_var_dict);
        optcontrol.use_sparse_solver = flag_sparse;
    }

    if (!fitting_var_dict["LASSO_DNORM"].empty()) {
        optcontrol.displacement_scaling_factor
            = boost::lexical_cast<double>(fitting_var_dict["LASSO_DNORM"]);
    }
    if (!fitting_var_dict["LASSO_ALPHA"].empty()) {
        optcontrol.l1_alpha = boost::lexical_cast<double>(fitting_var_dict["LASSO_ALPHA"]);
    }
    if (!fitting_var_dict["LASSO_MINALPHA"].empty()) {
        optcontrol.l1_alpha_min = boost::lexical_cast<double>(fitting_var_dict["LASSO_MINALPHA"]);
    }
    if (!fitting_var_dict["LASSO_MAXALPHA"].empty()) {
        optcontrol.l1_alpha_max = boost::lexical_cast<double>(fitting_var_dict["LASSO_MAXALPHA"]);
    }
    if (!fitting_var_dict["LASSO_NALPHA"].empty()) {
        optcontrol.num_l1_alpha = boost::lexical_cast<int>(fitting_var_dict["LASSO_NALPHA"]);
    }
    if (!fitting_var_dict["LASSO_TOL"].empty()) {
        optcontrol.tolerance_iteration = boost::lexical_cast<double>(fitting_var_dict["LASSO_TOL"]);
    }
    if (!fitting_var_dict["LASSO_MAXITER"].empty()) {
        optcontrol.maxnum_iteration = boost::lexical_cast<int>(fitting_var_dict["LASSO_MAXITER"]);
    }
    if (!fitting_var_dict["LASSO_CV"].empty()) {
        optcontrol.cross_validation_mode = boost::lexical_cast<int>(fitting_var_dict["LASSO_CV"]);
    }
    if (!fitting_var_dict["LASSO_CVSET"].empty()) {
        optcontrol.nset_cross_validation = boost::lexical_cast<int>(fitting_var_dict["LASSO_CVSET"]);
    }
    if (!fitting_var_dict["LASSO_FREQ"].empty()) {
        optcontrol.output_frequency = boost::lexical_cast<int>(fitting_var_dict["LASSO_FREQ"]);
    }
    if (!fitting_var_dict["STANDARDIZE"].empty()) {
        optcontrol.standardize = boost::lexical_cast<int>(fitting_var_dict["STANDARDIZE"]);
    }
    if (!fitting_var_dict["SOLUTION_PATH"].empty()) {
        optcontrol.save_solution_path = boost::lexical_cast<int>(fitting_var_dict["SOLUTION_PATH"]);
    }
    if (!fitting_var_dict["DEBIAS_OLS"].empty()) {
        optcontrol.debiase_after_l1opt = boost::lexical_cast<int>(fitting_var_dict["DEBIAS_OLS"]);
    }
    if (!fitting_var_dict["L1_RATIO"].empty()) {
        optcontrol.l1_ratio = boost::lexical_cast<double>(fitting_var_dict["L1_RATIO"]);
    }

    ndata_test = ndata;

    if (!fitting_var_dict["NDATA_TEST"].empty()) {
        ndata_test = boost::lexical_cast<int>(fitting_var_dict["NDATA_TEST"]);
    }

    if (fitting_var_dict["NSTART_TEST"].empty()) {
        nstart_test = 1;
    } else {
        nstart_test = boost::lexical_cast<int>(fitting_var_dict["NSTART_TEST"]);
    }
    if (fitting_var_dict["NEND_TEST"].empty()) {
        nend_test = 1;
    } else {
        nend_test = boost::lexical_cast<int>(fitting_var_dict["NEND_TEST"]);
    }

    if (fitting_var_dict["DFILE_TEST"].empty()) {
        dfile_test = dfile;
    } else {
        dfile_test = fitting_var_dict["DFILE_TEST"];
    }

    if (fitting_var_dict["FFILE_TEST"].empty()) {
        ffile_test = ffile;
    } else {
        ffile_test = fitting_var_dict["FFILE_TEST"];
    }

    if (ndata_test <= 0 || nstart_test <= 0 || nend_test <= 0
        || nstart_test > ndata_test || nend_test > ndata_test || nstart_test > nend_test) {
        exit("parse_fitting_vars",
             "ndata_test, nstart_test, nend_test are not consistent with each other");
    }

    input_setter->set_optimize_vars(alm,
                                    ndata,
                                    nstart,
                                    nend,
                                    skip_s,
                                    skip_e,
                                    dfile,
                                    ffile,
                                    ndata_test,
                                    nstart_test,
                                    nend_test,
                                    dfile_test,
                                    ffile_test,
                                    optcontrol);

    input_setter->set_constraint_vars(alm,
                                      constraint_flag,
                                      rotation_axis,
                                      fc2_file,
                                      fc3_file,
                                      fix_harmonic,
                                      fix_cubic);

    fitting_var_dict.clear();
}

void InputParser::parse_atomic_positions(ALM *alm)
{
    int i, j;
    std::string line, line_wo_comment;
    std::string str_tmp;
    std::string::size_type pos_first_comment_tag;
    std::vector<std::string> str_v, pos_line;
    double (*xeq)[3];
    int *kd;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }


    str_v.clear();

    if (from_stdin) {

        while (std::getline(std::cin, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\n\r "));
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_v.push_back(line_wo_comment);
        }

    } else {

        while (std::getline(ifs_input, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\n\r "));
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_v.push_back(line_wo_comment);
        }
    }


    if (str_v.size() != nat) {
        exit("parse_atomic_positions",
             "The number of entries for atomic positions should be NAT");
    }

    allocate(xeq, nat);
    allocate(kd, nat);


    for (i = 0; i < nat; ++i) {

        split_str_by_space(str_v[i], pos_line);

        if (pos_line.size() == 4) {
            try {
                kd[i] = boost::lexical_cast<int>(pos_line[0]);
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                exit("parse_atomic_positions",
                     "Invalid entry for the &position field at line ",
                     i + 1);
            }

            for (j = 0; j < 3; ++j) {
                xeq[i][j] = boost::lexical_cast<double>(pos_line[j + 1]);
            }

        } else {
            exit("parse_atomic_positions",
                 "Bad format for &position region");
        }
    }


    input_setter->set_atomic_positions(nat, kd, xeq);

    deallocate(xeq);
    deallocate(kd);
    pos_line.clear();
    str_v.clear();
}

void InputParser::parse_cutoff_radii()
{
    std::string line, line_wo_comment;
    std::string::size_type pos_first_comment_tag;
    std::vector<std::string> str_cutoff;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    str_cutoff.clear();

    if (from_stdin) {

        while (std::getline(std::cin, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_cutoff.push_back(line_wo_comment);
        }
    } else {

        while (std::getline(ifs_input, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_cutoff.push_back(line_wo_comment);
        }

    }

    int i, j, k;
    int order;
    std::vector<std::string> cutoff_line;
    std::set<std::string> element_allowed;
    std::vector<std::string> str_pair;
    std::map<std::string, int> kd_map;

    double cutoff_tmp;
    double ***cutoff_radii_tmp;
    bool ***undefined_cutoff;

    allocate(undefined_cutoff, maxorder, nkd, nkd);

    for (order = 0; order < maxorder; ++order) {
        for (i = 0; i < nkd; ++i) {
            for (j = 0; j < nkd; ++j) {
                undefined_cutoff[order][i][j] = true;
            }
        }
    }

    allocate(cutoff_radii_tmp, maxorder, nkd, nkd);

    element_allowed.clear();

    for (i = 0; i < nkd; ++i) {
        element_allowed.insert(kdname[i]);
        kd_map.insert(std::map<std::string, int>::value_type(kdname[i], i));
    }

    element_allowed.insert("*");
    kd_map.insert(std::map<std::string, int>::value_type("*", -1));

    for (std::vector<std::string>::const_iterator it = str_cutoff.begin();
         it != str_cutoff.end(); ++it) {

        split_str_by_space(*it, cutoff_line);

        if (cutoff_line.size() < maxorder + 1) {
            exit("parse_cutoff_radii",
                 "Invalid format for &cutoff entry");
        }

        boost::split(str_pair, cutoff_line[0], boost::is_any_of("-"));

        if (str_pair.size() != 2) {
            exit("parse_cutoff_radii2",
                 "Invalid format for &cutoff entry");
        }

        for (i = 0; i < 2; ++i) {
            if (element_allowed.find(str_pair[i]) == element_allowed.end()) {
                exit("parse_cutoff_radii2",
                     "Invalid format for &cutoff entry");
            }
        }

        const auto ikd = kd_map[str_pair[0]];
        const auto jkd = kd_map[str_pair[1]];

        for (order = 0; order < maxorder; ++order) {
            // Accept any strings starting with 'N' or 'n' as 'None'
            if ((cutoff_line[order + 1][0] == 'N') || (cutoff_line[order + 1][0] == 'n')) {
                // Minus value for cutoff radius.
                // This is a flag for neglecting cutoff radius
                cutoff_tmp = -1.0;
            } else {
                cutoff_tmp = boost::lexical_cast<double>(cutoff_line[order + 1]);
            }

            if (ikd == -1 && jkd == -1) {
                for (i = 0; i < nkd; ++i) {
                    for (j = 0; j < nkd; ++j) {
                        cutoff_radii_tmp[order][i][j] = cutoff_tmp;
                        undefined_cutoff[order][i][j] = false;
                    }
                }
            } else if (ikd == -1) {
                for (i = 0; i < nkd; ++i) {
                    cutoff_radii_tmp[order][i][jkd] = cutoff_tmp;
                    cutoff_radii_tmp[order][jkd][i] = cutoff_tmp;
                    undefined_cutoff[order][i][jkd] = false;
                    undefined_cutoff[order][jkd][i] = false;
                }
            } else if (jkd == -1) {
                for (j = 0; j < nkd; ++j) {
                    cutoff_radii_tmp[order][j][ikd] = cutoff_tmp;
                    cutoff_radii_tmp[order][ikd][j] = cutoff_tmp;
                    undefined_cutoff[order][j][ikd] = false;
                    undefined_cutoff[order][ikd][j] = false;
                }
            } else {
                cutoff_radii_tmp[order][ikd][jkd] = cutoff_tmp;
                cutoff_radii_tmp[order][jkd][ikd] = cutoff_tmp;
                undefined_cutoff[order][ikd][jkd] = false;
                undefined_cutoff[order][jkd][ikd] = false;
            }
        }
    }
    element_allowed.clear();
    str_cutoff.clear();

    for (order = 0; order < maxorder; ++order) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                if (undefined_cutoff[order][j][k]) {
                    std::cout << " Cutoff radius for " << std::setw(3)
                        << order + 2 << "th-order terms" << std::endl;
                    std::cout << " are not defined between elements " << std::setw(3) << j + 1
                        << " and " << std::setw(3) << k + 1 << std::endl;
                    exit("parse_cutoff_radii", "Incomplete cutoff radii");
                }
            }
        }
    }
    deallocate(undefined_cutoff);

    input_setter->set_cutoff_radii(maxorder,
                                   nkd,
                                   cutoff_radii_tmp);

    deallocate(cutoff_radii_tmp);
}

void InputParser::get_var_dict(const std::string keywords,
                               std::map<std::string, std::string> &var_dict)
{
    std::string line, key, val;
    std::string line_wo_comment;
    std::string::size_type pos_first_comment_tag;
    std::vector<std::string> str_entry, str_varval;

    std::set<std::string> keyword_set;

    boost::split(keyword_set, keywords, boost::is_space());

    var_dict.clear();

    if (from_stdin) {
        while (std::getline(std::cin, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            // Split the input line by ';'

            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));

            for (auto &it : str_entry) {

                // Split the input entry by '='

                auto str_tmp = boost::trim_copy(it);

                if (!str_tmp.empty()) {

                    boost::split(str_varval, str_tmp, boost::is_any_of("="));

                    if (str_varval.size() != 2) {
                        exit("get_var_dict", "Unacceptable format");
                    }

                    key = boost::to_upper_copy(boost::trim_copy(str_varval[0]));
                    val = boost::trim_copy(str_varval[1]);

                    if (keyword_set.find(key) == keyword_set.end()) {
                        std::cout << "Could not recognize the variable " << key << std::endl;
                        exit("get_var_dict", "Invalid variable found");
                    }

                    if (var_dict.find(key) != var_dict.end()) {
                        std::cout << "Variable " << key << " appears twice in the input file." << std::endl;
                        exit("get_var_dict", "Redundant input parameter");
                    }

                    // If everything is OK, add the variable and the corresponding value
                    // to the dictionary.

                    var_dict.insert(std::map<std::string, std::string>::value_type(key, val));
                }
            }
        }

    } else {

        while (std::getline(ifs_input, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            // Split the input line by ';'

            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));

            for (auto &it : str_entry) {

                // Split the input entry by '='

                std::string str_tmp = boost::trim_copy(it);

                if (!str_tmp.empty()) {

                    boost::split(str_varval, str_tmp, boost::is_any_of("="));

                    if (str_varval.size() != 2) {
                        exit("get_var_dict", "Unacceptable format");
                    }

                    key = boost::to_upper_copy(boost::trim_copy(str_varval[0]));
                    val = boost::trim_copy(str_varval[1]);

                    if (keyword_set.find(key) == keyword_set.end()) {
                        std::cout << "Could not recognize the variable "
                            << key << std::endl;
                        exit("get_var_dict", "Invalid variable found");
                    }

                    if (var_dict.find(key) != var_dict.end()) {
                        std::cout << "Variable " << key
                            << " appears twice in the input file." << std::endl;
                        exit("get_var_dict", "Redundant input parameter");
                    }

                    // If everything is OK, add the variable and the corresponding value
                    // to the dictionary.

                    var_dict.insert(std::map<std::string, std::string>::value_type(key, val));
                }
            }
        }
    }

    keyword_set.clear();
}


int InputParser::locate_tag(const std::string key)
{
    auto ret = 0;
    std::string line;

    if (from_stdin) {

        std::cin.clear();
        std::cin.seekg(0, std::ios_base::beg);

        while (std::cin >> line) {
            boost::to_lower(line);
            if (line == key) {
                ret = 1;
                break;
            }
        }
        return ret;

    } else {

        ifs_input.clear();
        ifs_input.seekg(0, std::ios_base::beg);

        while (ifs_input >> line) {
            boost::to_lower(line);
            if (line == key) {
                ret = 1;
                break;
            }
        }
        return ret;
    }
}

bool InputParser::is_endof_entry(const std::string str) const
{
    return str[0] == '/';
}

void InputParser::split_str_by_space(const std::string str,
                                     std::vector<std::string> &str_vec) const
{
    std::string str_tmp;
    std::istringstream is(str);

    str_vec.clear();

    while (1) {
        str_tmp.clear();
        is >> str_tmp;
        if (str_tmp.empty()) {
            break;
        }
        str_vec.push_back(str_tmp);
    }
    str_tmp.clear();
}

template <typename T>
void InputParser::assign_val(T &val,
                             const std::string key,
                             std::map<std::string, std::string> dict)
{
    // Assign a value to the variable "key" using the boost::lexica_cast.

    if (!dict[key].empty()) {
        try {
            val = boost::lexical_cast<T>(dict[key]);
        }
        catch (std::exception &e) {
            std::cout << e.what() << std::endl;
            std::string str_tmp = "Invalid entry for the " + key + " tag.\n";
            str_tmp += " Please check the input value.";
            exit("assign_val", str_tmp.c_str());
        }
    }
}
