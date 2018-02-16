/*
 constraint.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "constraint.h"
#include "combination.h"
#include "constants.h"
#include "error.h"
#include "fcs.h"
#include "interaction.h"
#include "mathfunctions.h"
#include "memory.h"
#include "rref.h"
#include "symmetry.h"
#include "system.h"
#include "timer.h"
#include "xml_parser.h"
#include <iostream>
#include <iomanip>
#include <boost/bimap.hpp>
#include <algorithm>
#include <unordered_set>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

Constraint::Constraint()
{
    set_default_variables();
}

Constraint::~Constraint()
{
    deallocate_variables();
}

void Constraint::set_default_variables()
{
    constraint_mode = 1;
    rotation_axis = "";
    fix_harmonic = false;
    fix_cubic = false;
    constraint_algebraic = 0;
    fc2_file = "";
    fc3_file = "";
    exist_constraint = false;
    extra_constraint_from_symmetry = false;
    const_mat = nullptr;
    const_rhs = nullptr;
    const_symmetry = nullptr;
    const_fix = nullptr;
    const_relate = nullptr;
    const_relate_rotation = nullptr;
    index_bimap = nullptr;
    number_of_constraints = 0;
    tolerance_constraint = eps6;
}

void Constraint::deallocate_variables()
{
    if (const_symmetry) {
        deallocate(const_symmetry);
    }
    if (const_fix) {
        deallocate(const_fix);
    }
    if (const_relate) {
        deallocate(const_relate);
    }
    if (const_relate_rotation) {
        deallocate(const_relate_rotation);
    }
    if (index_bimap) {
        deallocate(index_bimap);
    }
    if (const_mat) {
        deallocate(const_mat);
    }
    if (const_rhs) {
        deallocate(const_rhs);
    }
}

void Constraint::setup(ALM *alm)
{
    alm->timer->start_clock("constraint");
    std::cout << " CONSTRAINT" << std::endl;
    std::cout << " ==========" << std::endl << std::endl;

    constraint_algebraic = constraint_mode / 10;
    constraint_mode = constraint_mode % 10;
    int maxorder = alm->interaction->maxorder;

    switch (constraint_mode) {
    case 0: // do nothing
        impose_inv_T = false;
        impose_inv_R = false;
        std::cout << "  ICONST = 0: Constraint for translational/rotational invariance" << std::endl;
        std::cout << "              will NOT be considered." << std::endl;
        break;
    case 1:
        impose_inv_T = true;
        impose_inv_R = false;
        std::cout << "  ICONST = 1: Constraints for translational invariance" << std::endl;
        std::cout << "              will be considered." << std::endl;
        break;
    case 2:
        impose_inv_T = true;
        impose_inv_R = true;
        exclude_last_R = true;
        std::cout << "  ICONST = 2: Constraints for translational and rotational invariance" << std::endl;
        std::cout << "              will be considered. Axis of rotation is " << rotation_axis << std::endl;
        std::cout << "              Rotational invariance of the maximum order will be neglected" << std::endl;
        break;
    case 3:
        impose_inv_T = true;
        impose_inv_R = true;
        exclude_last_R = false;
        std::cout << "  ICONST = 3: Constraints for translational and rotational invariance" << std::endl;
        std::cout << "              will be considered. Axis of rotation is " << rotation_axis << std::endl;
        break;
    default:
        exit("Constraint::setup", "invalid constraint_mode", constraint_mode);
        break;
    }

    std::cout << std::endl;


    if (const_symmetry) {
        deallocate(const_symmetry);
    }
    if (const_fix) {
        deallocate(const_fix);
    }
    if (const_relate) {
        deallocate(const_relate);
    }
    if (index_bimap) {
        deallocate(index_bimap);
    }

    allocate(const_fix, maxorder);
    allocate(const_relate, maxorder);
    allocate(index_bimap, maxorder);
    allocate(const_symmetry, maxorder);
    allocate(const_translation, maxorder);
    allocate(const_rotation_self, maxorder);
    allocate(const_rotation_cross, maxorder);
    allocate(const_self, maxorder);

    for (int order = 0; order < maxorder; ++order) {
        const_translation[order].clear();
        const_rotation_self[order].clear();
        const_rotation_cross[order].clear();
        const_self[order].clear();
        const_fix[order].clear();
    }

    if (fix_harmonic) {

        std::cout << "  FC2XML is given : Harmonic force constants will be " << std::endl;
        std::cout << "                    fixed to the values given in " << fc2_file << std::endl;
        std::cout << std::endl;

        fix_forceconstants_to_file(0,
                                   alm->symmetry,
                                   alm->fcs,
                                   fc2_file,
                                   const_fix[0]);
    }

    fix_cubic = fix_cubic & (alm->interaction->maxorder > 1);
    if (fix_cubic) {

        std::cout << "  FC3XML is given : Cubic force constants will be " << std::endl;
        std::cout << "                    fixed to the values given in " << fc3_file << std::endl;
        std::cout << std::endl;

        fix_forceconstants_to_file(1,
                                   alm->symmetry,
                                   alm->fcs,
                                   fc3_file,
                                   const_fix[1]);
    }

    generate_symmetry_constraint_in_cartesian(alm->system->supercell.number_of_atoms,
                                              alm->symmetry,
                                              alm->interaction,
                                              alm->fcs,
                                              const_symmetry);

    extra_constraint_from_symmetry = false;

    for (int order = 0; order < alm->interaction->maxorder; ++order) {
        if (!const_symmetry[order].empty()) extra_constraint_from_symmetry = true;
    }

    if (impose_inv_T) {
        generate_translational_constraint(alm->system->supercell,
                                          alm->symmetry,
                                          alm->interaction,
                                          alm->fcs,
                                          const_translation);
    }

    if (impose_inv_R) {
        generate_rotational_constraint(alm->system,
                                       alm->symmetry,
                                       alm->interaction,
                                       alm->fcs,
                                       const_rotation_self,
                                       const_rotation_cross);
    }

    // double *arr_tmp;

    // Merge intra-order constrants and do reduction 

    for (int order = 0; order < maxorder; ++order) {

        int nparam = alm->fcs->nequiv[order].size();

        const_self[order].reserve(
            const_translation[order].size()
            + const_rotation_self[order].size()
            + const_symmetry[order].size());

        const_self[order].insert(const_self[order].end(),
                                 const_translation[order].begin(),
                                 const_translation[order].end());

        const_self[order].insert(const_self[order].end(),
                                 const_rotation_self[order].begin(),
                                 const_rotation_self[order].end());

        const_self[order].insert(const_self[order].end(),
                                 const_symmetry[order].begin(),
                                 const_symmetry[order].end());

        remove_redundant_rows(nparam, const_self[order], eps8);

        //allocate(arr_tmp, nparam);

        //for (const auto &e : const_translation[order]) {
        //    for (int i = 0; i < nparam; ++i) {
        //        arr_tmp[i] = e.w_const[i];
        //    }
        //    const_self[order].emplace_back(nparam, arr_tmp);
        //}

        //if (!const_rotation_self[order].empty()) {
        //    for (const auto &e : const_rotation_self[order]) {
        //        for (int i = 0; i < nparam; ++i) {
        //            arr_tmp[i] = e.w_const[i];
        //        }
        //        const_self[order].emplace_back(nparam, arr_tmp);
        //    }
        //    remove_redundant_rows(nparam, const_self[order], eps8);
        //}

        //if (!const_symmetry[order].empty()) {
        //    for (const auto &e : const_symmetry[order]) {
        //        for (int i = 0; i < nparam; ++i) {
        //            arr_tmp[i] = e.w_const[i];
        //        }
        //        const_self[order].emplace_back(nparam, arr_tmp);
        //    }
        //    remove_redundant_rows(nparam, const_self[order], eps8);
        //}

        //deallocate(arr_tmp);
    }

    if (constraint_algebraic) {

        get_mapping_constraint(maxorder,
                               alm->fcs->nequiv,
                               const_self,
                               const_fix,
                               const_relate,
                               index_bimap);


    } else {

        int Pmax = 0;
        int nparams = 0;
        for (int order = 0; order < maxorder; ++order) {
            Pmax += const_self[order].size()
                + const_rotation_cross[order].size();
        }
        if (fix_harmonic) {
            Pmax -= const_self[0].size();
            Pmax += alm->fcs->nequiv[0].size();
        }
        if (fix_cubic) {
            Pmax -= const_self[1].size();
            Pmax += alm->fcs->nequiv[1].size();
        }
        for (int order = 0; order < maxorder; ++order) {
            nparams += alm->fcs->nequiv[order].size();
        }

        if (const_mat) {
            deallocate(const_mat);
        }
        allocate(const_mat, Pmax, nparams);

        if (const_rhs) {
            deallocate(const_rhs);
        }
        allocate(const_rhs, Pmax);

        calc_constraint_matrix(maxorder,
                               alm->fcs->nequiv,
                               nparams,
                               number_of_constraints);
    }

    exist_constraint
        = impose_inv_T
        || fix_harmonic
        || fix_cubic
        || extra_constraint_from_symmetry;

    if (exist_constraint) {

        int order;

        if (impose_inv_T || impose_inv_R) {
            std::cout << "  Number of constraints [T-inv, R-inv (self), R-inv (cross)]:" << std::endl;
            for (order = 0; order < maxorder; ++order) {
                std::cout << "   " << std::setw(8) << alm->interaction->str_order[order];
                std::cout << std::setw(5) << const_translation[order].size();
                std::cout << std::setw(5) << const_rotation_self[order].size();
                std::cout << std::setw(5) << const_rotation_cross[order].size();
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        if (extra_constraint_from_symmetry) {
            std::cout << "  There are constraints from crystal symmetry." << std::endl;
            std::cout << "  The number of such constraints for each order:" << std::endl;
            for (order = 0; order < maxorder; ++order) {
                std::cout << "   " << std::setw(8) << alm->interaction->str_order[order];
                std::cout << std::setw(5) << const_symmetry[order].size();
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        if (extra_constraint_from_symmetry) {
            std::cout << "  Constraints of T-inv, R-inv (self), and those from crystal symmetry are merged."
                << std::endl;
        } else {
            std::cout << "  Constraints of T-inv and R-inv (self) are merged." << std::endl;
        }
        std::cout << "  If there are redundant constraints, they are removed in this process." << std::endl;
        std::cout << std::endl;
        std::cout << "  Number of inequivalent constraints (self, cross) : " << std::endl;

        for (order = 0; order < maxorder; ++order) {
            std::cout << "   " << std::setw(8) << alm->interaction->str_order[order];
            std::cout << std::setw(5) << const_self[order].size();
            std::cout << std::setw(5) << const_rotation_cross[order].size();
            std::cout << std::endl;
        }
        std::cout << std::endl;

        if (constraint_algebraic) {

            std::cout << "  ICONST >= 10 : Constraints will be considered algebraically."
                << std::endl << std::endl;

            if (impose_inv_R) {
                std::cout << "  WARNING : Inter-order constraints for rotational invariance will be neglected."
                    << std::endl;
            }

            for (order = 0; order < maxorder; ++order) {
                std::cout << "  Number of free" << std::setw(9) << alm->interaction->str_order[order]
                    << " FCs : " << index_bimap[order].size() << std::endl;
            }
            std::cout << std::endl;

        } else {

            std::cout << "  Total number of constraints = " << number_of_constraints << std::endl << std::endl;

        }
    }

    deallocate(const_translation);
    const_translation = nullptr;
    deallocate(const_rotation_self);
    const_rotation_self = nullptr;
    deallocate(const_rotation_cross);
    const_rotation_cross = nullptr;
    deallocate(const_self);
    const_self = nullptr;

    alm->timer->print_elapsed();
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    alm->timer->stop_clock("constraint");
}

void Constraint::calc_constraint_matrix(const int maxorder,
                                        std::vector<int> *nequiv,
                                        const int nparams,
                                        int &nconst)
{
    int i, j;
    int order;
    double *arr_tmp;
    std::vector<ConstraintClass> const_total;

    const_total.clear();
    allocate(arr_tmp, nparams);

    int nshift = 0;

    for (order = 0; order < maxorder; ++order) {
        int nelems = nequiv[order].size();

        if (const_fix[order].empty()) {

            for (i = 0; i < nelems; ++i) arr_tmp[i] = 0.0;

            for (auto &p : const_self[order]) {
                for (i = 0; i < nelems; ++i) {
                    arr_tmp[nshift + i] = p.w_const[i];
                }
                const_total.emplace_back(nparams, arr_tmp);
            }
        }
        nshift += nelems;
    }

    auto nconst1 = const_total.size();

    // Inter-order constraints
    int nshift2 = 0;
    for (order = 0; order < maxorder; ++order) {
        if (order > 0) {
            int nparam2 = nequiv[order - 1].size() + nequiv[order].size();

            if (const_fix[order - 1].empty() && const_fix[order].empty()) {

                for (i = 0; i < nparams; ++i) arr_tmp[i] = 0.0;

                for (auto &p : const_rotation_cross[order]) {
                    for (i = 0; i < nparam2; ++i) {
                        arr_tmp[nshift2 + i] = p.w_const[i];
                    }
                    const_total.emplace_back(nparams, arr_tmp);
                }
            }

            nshift2 += nequiv[order - 1].size();
        }
    }
    deallocate(arr_tmp);

    if (nconst1 != const_total.size())
        remove_redundant_rows(nparams, const_total, eps8);

    nconst = const_total.size();

    if (fix_harmonic) nconst += nequiv[0].size();
    if (fix_cubic) nconst += nequiv[1].size();

    /*  if (fix_harmonic) {
          std::cout << "  Harmonic force constants will be fixed to the values " << std::endl;
          std::cout << "  of the reference " << fc2_file << std::endl;
          std::cout << "  Constraint for HARMONIC IFCs will be updated accordingly."
              << std::endl << std::endl;
          nconst += nequiv[0].size();
      }
  
      if (fix_cubic) {
          std::cout << "  Cubic force constants will be fixed to the values " << std::endl;
          std::cout << "  of the reference " << fc3_file << std::endl;
          std::cout << "  Constraint for ANHARM3 IFCs will be updated accordingly."
              << std::endl << std::endl;
          nconst += nequiv[1].size();
      }*/

    for (i = 0; i < nconst; ++i) {
        for (j = 0; j < nparams; ++j) {
            const_mat[i][j] = 0.0;
        }
        const_rhs[i] = 0.0;
    }

    int irow = 0;
    int icol = 0;

    int ishift = 0;


    if (fix_harmonic) {

        for (const auto &p : const_fix[0]) {
            i = p.p_index_target;
            const_mat[i][i] = 1.0;
            const_rhs[i] = p.val_to_fix;
        }

        irow += const_fix[0].size();
        icol += const_fix[0].size();;
        ishift += const_fix[0].size();
    }

    if (fix_cubic && maxorder > 1) {

        int ishift2 = nequiv[0].size();

        for (const auto &p : const_fix[1]) {
            i = p.p_index_target;
            const_mat[i + ishift][i + ishift2] = 1.0;
            const_rhs[i + ishift] = p.val_to_fix;
        }

        irow += const_fix[1].size();
        icol += const_fix[1].size();
    }

    for (auto &p : const_total) {
        for (i = 0; i < nparams; ++i) {
            const_mat[irow][i] = p.w_const[i];
        }
        ++irow;
    }
    const_total.clear();
}


void Constraint::get_mapping_constraint(const int nmax,
                                        std::vector<int> *nequiv,
                                        std::vector<ConstraintClass> *const_in,
                                        std::vector<ConstraintTypeFix> *const_fix_out,
                                        std::vector<ConstraintTypeRelate> *const_relate_out,
                                        boost::bimap<int, int> *index_bimap_out)
{
    // If const_fix_out[order] is not empty as input, it assumes that fix_forceconstant[order] is true.
    // In this case, const_fix_out[order] is not updated.

    int order;
    unsigned int i;

    int nparam;
    for (order = 0; order < nmax; ++order) {

        nparam = nequiv[order].size();

        if (const_fix_out[order].empty()) {

            int p_index_target;
            std::vector<double> alpha_tmp;
            std::vector<unsigned int> p_index_tmp;

            for (auto p = const_in[order].rbegin(); p != const_in[order].rend(); ++p) {
                p_index_target = -1;
                for (i = 0; i < nparam; ++i) {
                    if (std::abs((*p).w_const[i]) > tolerance_constraint) {
                        p_index_target = i;
                        break;
                    }
                }

                if (p_index_target == -1) {
                    exit("get_mapping_constraint",
                         "No finite entry found in the constraint.");
                }

                alpha_tmp.clear();
                p_index_tmp.clear();

                for (i = p_index_target + 1; i < nparam; ++i) {
                    if (std::abs((*p).w_const[i]) > tolerance_constraint) {
                        alpha_tmp.push_back((*p).w_const[i]);
                        p_index_tmp.push_back(i);
                    }
                }


                if (!alpha_tmp.empty()) {
                    const_relate_out[order].emplace_back(p_index_target,
                                                         alpha_tmp, p_index_tmp);
                } else {
                    const_fix_out[order].emplace_back(p_index_target, 0.0);
                }
            }
        }
    }

    std::vector<int> *has_constraint;
    allocate(has_constraint, nmax);

    for (order = 0; order < nmax; ++order) {

        nparam = nequiv[order].size();

        for (i = 0; i < nparam; ++i) {
            has_constraint[order].push_back(0);
        }

        for (i = 0; i < const_fix_out[order].size(); ++i) {
            has_constraint[order][const_fix_out[order][i].p_index_target] = 1;
        }

        for (i = 0; i < const_relate_out[order].size(); ++i) {
            has_constraint[order][const_relate_out[order][i].p_index_target] = 2;
        }
    }

    int icount;

    for (order = 0; order < nmax; ++order) {
        nparam = nequiv[order].size();

        icount = 0;
        for (i = 0; i < nparam; ++i) {
            if (!has_constraint[order][i]) {
                index_bimap_out[order].insert(
                    boost::bimap<int, int>::value_type(icount, i));
                ++icount;
            }
        }
    }

    deallocate(has_constraint);
}

void Constraint::generate_symmetry_constraint_in_cartesian(const int nat,
                                                           Symmetry *symmetry,
                                                           Interaction *interaction,
                                                           Fcs *fcs,
                                                           std::vector<ConstraintClass> *const_out)

{
    // Create constraint matrices arising from the crystal symmetry.

    int maxorder = interaction->maxorder;
    bool has_constraint_from_symm = false;
    std::vector<std::vector<double>> const_tmp;

    for (int isym = 0; isym < symmetry->nsym; ++isym) {
        if (!symmetry->SymmData[isym].compatible_with_cartesian) {
            has_constraint_from_symm = true;
            break;
        }
    }

    if (has_constraint_from_symm) {
        std::cout << "  Generating constraints from crystal symmetry ..." << std::endl;
    }

    for (int order = 0; order < maxorder; ++order) {
        if (has_constraint_from_symm) {
            std::cout << "   " << std::setw(8) << interaction->str_order[order] << " ...";
        }
        fcs->get_constraint_symmetry(nat,
                                     symmetry,
                                     order,
                                     interaction->cluster_list[order],
                                     "Cartesian",
                                     fcs->fc_table[order],
                                     fcs->nequiv[order].size(),
                                     tolerance_constraint,
                                     const_tmp);

        for (auto &it : const_tmp) {
            const_out[order].emplace_back(ConstraintClass(it));
        }

        if (has_constraint_from_symm) {
            std::cout << " done." << std::endl;
        }
    }
    if (has_constraint_from_symm) {
        std::cout << "  Finished !" << std::endl << std::endl;
    }
}


void Constraint::generate_translational_constraint(const Cell &supercell,
                                                   Symmetry *symmetry,
                                                   Interaction *interaction,
                                                   Fcs *fcs,
                                                   std::vector<ConstraintClass> *const_out)
{
    // Create constraint matrix for the translational invariance (aka acoustic sum rule).

    std::cout << "  Generating constraints for translational invariance ..." << std::endl;

    for (int order = 0; order < interaction->maxorder; ++order) {

        std::cout << "   " << std::setw(8) << interaction->str_order[order] << " ...";

        int nparams = fcs->nequiv[order].size();

        if (nparams == 0) {
            std::cout << "  No parameters! Skipped." << std::endl;
            continue;
        }

        get_constraint_translation(supercell,
                                   symmetry,
                                   interaction,
                                   fcs,
                                   order,
                                   fcs->fc_table[order],
                                   fcs->nequiv[order].size(),
                                   const_out[order]);


        std::cout << " done." << std::endl;
    }

    std::cout << "  Finished !" << std::endl << std::endl;
}


void Constraint::get_constraint_translation(const Cell &supercell,
                                            Symmetry *symmetry,
                                            Interaction *interaction,
                                            Fcs *fcs,
                                            const int order,
                                            const std::vector<FcProperty> &fc_table,
                                            const int nparams,
                                            std::vector<ConstraintClass> &const_out)
{
    // Generate equality constraint for the acoustic sum rule.

    int i, j;
    int iat, jat, icrd, jcrd;
    int idata;
    int loc_nonzero;

    int *ind;
    int *intarr, *intarr_copy;
    int **xyzcomponent;

    int ixyz, nxyz;
    int natmin = symmetry->nat_prim;
    int nat = supercell.number_of_atoms;

    unsigned int isize;
    double *arr_constraint;

    std::vector<int> data;
    std::unordered_set<FcProperty> list_found;
    std::unordered_set<FcProperty>::iterator iter_found;
    std::vector<std::vector<int>> data_vec;
    std::vector<FcProperty> list_vec;
    std::vector<int> const_now;
    std::vector<std::vector<int>> const_mat;


    if (order < 0) return;

    const_mat.clear();

    if (nparams == 0) return;

    allocate(ind, order + 2);

    // Create force constant table for search

    list_found.clear();

    for (const auto &p : fc_table) {
        for (i = 0; i < order + 2; ++i) {
            ind[i] = p.elems[i];
        }
        if (list_found.find(FcProperty(order + 2, p.sign,
                                       ind, p.mother)) != list_found.end()) {
            exit("get_constraint_translation", "Duplicate interaction list found");
        }
        list_found.insert(FcProperty(order + 2, p.sign,
                                     ind, p.mother));
    }

    deallocate(ind);

    // Generate xyz component for each order

    nxyz = static_cast<int>(std::pow(static_cast<double>(3), order + 1));
    allocate(xyzcomponent, nxyz, order + 1);
    fcs->get_xyzcomponent(order + 1, xyzcomponent);

    allocate(intarr, order + 2);
    allocate(intarr_copy, order + 2);

    const_now.resize(nparams);

    for (i = 0; i < natmin; ++i) {

        iat = symmetry->map_p2s[i][0];

        // Generate atom pairs for each order

        if (order == 0) {
            for (icrd = 0; icrd < 3; ++icrd) {

                intarr[0] = 3 * iat + icrd;

                for (jcrd = 0; jcrd < 3; ++jcrd) {

                    // Reset the temporary array for another constraint
                    for (j = 0; j < nparams; ++j) const_now[j] = 0;

                    for (jat = 0; jat < 3 * nat; jat += 3) {
                        intarr[1] = jat + jcrd;

                        iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                intarr, 1));

                        //  If found an IFC
                        if (iter_found != list_found.end()) {
                            // Round the coefficient to integer
                            const_now[(*iter_found).mother] += nint((*iter_found).sign);
                        }

                    }
                    // Add to the constraint list
                    if (!is_allzero(const_now, loc_nonzero)) {
                        if (const_now[loc_nonzero] < 0) {
                            for (j = 0; j < nparams; ++j) const_now[j] *= -1;
                        }
                        const_mat.push_back(const_now);
                    }
                }
            }

        } else {

            // Anharmonic cases

            std::vector<int> intlist(interaction->interaction_pair[order][i]);
            std::sort(intlist.begin(), intlist.end());

            data_vec.clear();
            // Generate data_vec that contains possible interaction clusters.
            // Each cluster contains (order + 1) atoms, and the last atom index
            // will be treated seperately below.
            CombinationWithRepetition<int> g2(intlist.begin(), intlist.end(), order);
            do {
                data = g2.now();

                intarr[0] = iat;

                for (isize = 0; isize < data.size(); ++isize) {
                    intarr[isize + 1] = data[isize];
                }

                if (interaction->satisfy_nbody_rule(order + 1, intarr, order)) {
                    if (interaction->is_incutoff(order + 1, intarr, order, supercell.kind)) {
                        // Add to list if the atoms interact with each other.
                        data_vec.push_back(data);
                    }
                }

            } while (g2.next());

            int ndata = data_vec.size();

            // Use openmp for acceleration if possible
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                int *intarr_omp, *intarr_copy_omp;

                allocate(intarr_omp, order + 2);
                allocate(intarr_copy_omp, order + 2);

                std::vector<std::vector<int>> const_omp;
                std::vector<int> data_omp;
                std::vector<int> const_now_omp;

                const_omp.clear();
                const_now_omp.resize(nparams);
#ifdef _OPENMP
#pragma omp for private(isize, ixyz, jcrd, j, jat, iter_found, loc_nonzero), schedule(guided), nowait
#endif
                for (idata = 0; idata < ndata; ++idata) {

                    data_omp = data_vec[idata];

                    intarr_omp[0] = iat;
                    for (isize = 0; isize < data_omp.size(); ++isize) {
                        intarr_omp[isize + 1] = data_omp[isize];
                    }

                    // Loop for xyz component
                    for (ixyz = 0; ixyz < nxyz; ++ixyz) {
                        // Loop for the xyz index of the last atom
                        for (jcrd = 0; jcrd < 3; ++jcrd) {

                            // Reset the temporary array for another constraint
                            for (j = 0; j < nparams; ++j) const_now_omp[j] = 0;

                            // Loop for the last atom index
                            for (jat = 0; jat < 3 * nat; jat += 3) {
                                intarr_omp[order + 1] = jat / 3;

                                if (interaction->satisfy_nbody_rule(order + 2, intarr_omp, order)) {
                                    for (j = 0; j < order + 1; ++j) {
                                        intarr_copy_omp[j] = 3 * intarr_omp[j] + xyzcomponent[ixyz][j];
                                    }
                                    intarr_copy_omp[order + 1] = jat + jcrd;

                                    sort_tail(order + 2, intarr_copy_omp);

                                    iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                            intarr_copy_omp, 1));
                                    if (iter_found != list_found.end()) {
                                        const_now_omp[(*iter_found).mother] += nint((*iter_found).sign);
                                    }

                                }
                            } // close loop jat

                            // Add the constraint to the private array
                            if (!is_allzero(const_now_omp, loc_nonzero)) {
                                if (const_now_omp[loc_nonzero] < 0) {
                                    for (j = 0; j < nparams; ++j) const_now_omp[j] *= -1;
                                }
                                const_omp.push_back(const_now_omp);
                            }
                        }
                    }
                    // sort-->uniq the array
                    std::sort(const_omp.begin(), const_omp.end());
                    const_omp.erase(std::unique(const_omp.begin(), const_omp.end()),
                                    const_omp.end());

                    // Merge vectors
#pragma omp critical
                    {
                        for (auto &it : const_omp) {
                            const_mat.push_back(it);
                        }
                    }
                    const_omp.clear();

                } // close idata (openmp main loop)

                deallocate(intarr_omp);
                deallocate(intarr_copy_omp);
            } // close openmp

            intlist.clear();
        } // close if

        // sort--> uniq the array (to save memory consumption)
        std::sort(const_mat.begin(), const_mat.end());
        const_mat.erase(std::unique(const_mat.begin(), const_mat.end()),
                        const_mat.end());
    } // close loop i

    deallocate(xyzcomponent);
    deallocate(intarr);
    deallocate(intarr_copy);

    // Copy to constraint class 

    const_out.clear();
    allocate(arr_constraint, nparams);
    for (auto it = const_mat.rbegin(); it != const_mat.rend(); ++it) {
        for (i = 0; i < (*it).size(); ++i) {
            arr_constraint[i] = static_cast<double>((*it)[i]);
        }
        const_out.emplace_back(ConstraintClass(nparams,
                                               arr_constraint));
    }
    deallocate(arr_constraint);
    const_mat.clear();

    // Transform the matrix into the reduced row echelon form
    remove_redundant_rows(nparams, const_out, eps8);
}


void Constraint::generate_rotational_constraint(System *system,
                                                Symmetry *symmetry,
                                                Interaction *interaction,
                                                Fcs *fcs,
                                                std::vector<ConstraintClass> *const_rotation_self,
                                                std::vector<ConstraintClass> *const_rotation_cross)
{
    // Create constraints for the rotational invariance

    std::cout << "  Generating constraints for rotational invariance ..." << std::endl;

#ifdef _DEBUG
    std::ofstream ofs_constraint;
    ofs_constraint.open("CONSTRAINT", std::ios::out);
#endif

    int i, j;
    int iat, jat;
    int icrd, jcrd;
    int order;
    int maxorder = interaction->maxorder;
    int natmin = symmetry->nat_prim;
    int mu, nu;
    int ixyz, nxyz, nxyz2;
    int mu_lambda, lambda;
    int levi_factor;

    int *ind;
    int **xyzcomponent, **xyzcomponent2;
    int *nparams, nparam_sub;
    int *interaction_index, *interaction_atom;
    int *interaction_tmp;

    double *arr_constraint;
    double *arr_constraint_self;

    bool valid_rotation_axis[3][3];

    double vec_for_rot[3];

    std::vector<int> interaction_list;

    std::unordered_set<FcProperty> list_found;
    std::unordered_set<FcProperty> list_found_last;
    std::unordered_set<FcProperty>::iterator iter_found;

    CombinationWithRepetition<int> g;

    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<InteractionCluster>::iterator iter_cluster;

    setup_rotation_axis(valid_rotation_axis);

    allocate(ind, maxorder + 1);
    allocate(nparams, maxorder);

    for (order = 0; order < maxorder; ++order) {

        nparams[order] = fcs->nequiv[order].size();

        if (order == 0) {
            std::cout << "   Constraints between " << std::setw(8)
                << "1st-order IFCs (which are zero) and "
                << std::setw(8) << interaction->str_order[order] << " ...";
            nparam_sub = nparams[order];
        } else {
            std::cout << "   Constraints between " << std::setw(8)
                << interaction->str_order[order - 1] << " and "
                << std::setw(8) << interaction->str_order[order] << " ...";
            nparam_sub = nparams[order] + nparams[order - 1];
        }

        allocate(arr_constraint, nparam_sub);
        allocate(arr_constraint_self, nparams[order]);
        allocate(interaction_atom, order + 2);
        allocate(interaction_index, order + 2);
        allocate(interaction_tmp, order + 2);

        if (order > 0) {
            list_found_last = list_found;
            nxyz = static_cast<int>(pow(static_cast<double>(3), order));
            allocate(xyzcomponent, nxyz, order);
            fcs->get_xyzcomponent(order, xyzcomponent);
        }

        list_found.clear();

        for (auto p = fcs->fc_table[order].begin(); p != fcs->fc_table[order].end(); ++p) {
            for (i = 0; i < order + 2; ++i) {
                ind[i] = (*p).elems[i];
            }
            list_found.insert(FcProperty(order + 2, (*p).sign,
                                         ind, (*p).mother));
        }

        for (i = 0; i < natmin; ++i) {

            iat = symmetry->map_p2s[i][0];

            interaction_atom[0] = iat;

            if (order == 0) {

                std::vector<int> interaction_list_now(interaction->interaction_pair[order][i]);
                std::sort(interaction_list_now.begin(), interaction_list_now.end());

                // Special treatment for harmonic force constants

                for (icrd = 0; icrd < 3; ++icrd) {

                    interaction_index[0] = 3 * iat + icrd;

                    for (mu = 0; mu < 3; ++mu) {

                        for (nu = 0; nu < 3; ++nu) {

                            if (!valid_rotation_axis[mu][nu]) continue;

                            // Clear history

                            for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;

                            for (auto &iter_list : interaction_list_now) {

                                jat = iter_list;
                                interaction_index[1] = 3 * jat + mu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                        interaction_index, 1));

                                atom_tmp.clear();
                                atom_tmp.push_back(jat);
                                cell_dummy.clear();
                                iter_cluster = interaction->interaction_cluster[order][i].find(
                                    InteractionCluster(atom_tmp, cell_dummy));

                                if (iter_cluster == interaction->interaction_cluster[order][i].end()) {
                                    exit("generate_rotational_constraint",
                                         "interaction not found ...");
                                } else {
                                    for (j = 0; j < 3; ++j) vec_for_rot[j] = 0.0;

                                    int nsize_equiv = (*iter_cluster).cell.size();

                                    for (j = 0; j < nsize_equiv; ++j) {
                                        for (int k = 0; k < 3; ++k) {
                                            vec_for_rot[k]
                                                += system->x_image[(*iter_cluster).cell[j][0]][jat][k];
                                        }
                                    }

                                    for (j = 0; j < 3; ++j) {
                                        vec_for_rot[j] /= static_cast<double>(nsize_equiv);
                                    }
                                }


                                if (iter_found != list_found.end()) {
                                    arr_constraint[(*iter_found).mother] += (*iter_found).sign * vec_for_rot[nu];
                                }

                                // Exchange mu <--> nu and repeat again. 
                                // Note that the sign is inverted (+ --> -) in the summation

                                interaction_index[1] = 3 * jat + nu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                        interaction_index, 1));
                                if (iter_found != list_found.end()) {
                                    arr_constraint[(*iter_found).mother]
                                        -= (*iter_found).sign * vec_for_rot[mu];
                                }
                            }

                            if (!is_allzero(nparam_sub, arr_constraint)) {
                                // Add to constraint list
                                const_rotation_self[order].emplace_back(
                                    ConstraintClass(nparam_sub, arr_constraint));
                            }

                        } // nu
                    } // mu
                }
            } else {

                // Constraint between different orders

                std::vector<int> interaction_list_now(interaction->interaction_pair[order][i]);
                std::vector<int> interaction_list_old(interaction->interaction_pair[order - 1][i]);
                std::sort(interaction_list_now.begin(), interaction_list_now.end());
                std::sort(interaction_list_old.begin(), interaction_list_old.end());

                for (icrd = 0; icrd < 3; ++icrd) {

                    interaction_index[0] = 3 * iat + icrd;

                    CombinationWithRepetition<int> g_now(interaction_list_now.begin(),
                                                         interaction_list_now.end(), order);
                    CombinationWithRepetition<int> g_old(interaction_list_old.begin(),
                                                         interaction_list_old.end(), order);

                    // m    -th order --> (m-1)-th order
                    // (m-1)-th order -->     m-th order
                    // 2-different directions to find all constraints

                    for (unsigned int direction = 0; direction < 2; ++direction) {

                        if (direction == 0) {
                            g = g_now;
                            interaction_list = interaction_list_now;
                        } else {
                            g = g_old;
                            interaction_list = interaction_list_old;
                        }

                        // Loop for the interacting pairs

                        do {
                            std::vector<int> data = g.now();

                            for (int idata = 0; idata < data.size(); ++idata)
                                interaction_atom[idata + 1] = data[idata];

                            for (ixyz = 0; ixyz < nxyz; ++ixyz) {

                                for (j = 0; j < order; ++j)
                                    interaction_index[j + 1] = 3 * interaction_atom[j + 1] + xyzcomponent[ixyz][j];

                                for (mu = 0; mu < 3; ++mu) {

                                    for (nu = 0; nu < 3; ++nu) {

                                        if (!valid_rotation_axis[mu][nu]) continue;

                                        // Search for a new constraint below

                                        for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;

                                        // Loop for m_{N+1}, a_{N+1}
                                        for (auto &iter_list : interaction_list) {
                                            jat = iter_list;

                                            interaction_atom[order + 1] = jat;
                                            if (!interaction->is_incutoff(order + 2,
                                                                          interaction_atom,
                                                                          order,
                                                                          system->supercell.kind))
                                                continue;

                                            atom_tmp.clear();

                                            for (j = 1; j < order + 2; ++j) {
                                                atom_tmp.push_back(interaction_atom[j]);
                                            }
                                            std::sort(atom_tmp.begin(), atom_tmp.end());

                                            iter_cluster = interaction->interaction_cluster[order][i].find(
                                                InteractionCluster(atom_tmp,
                                                                   cell_dummy));
                                            if (iter_cluster != interaction->interaction_cluster[order][i].end()) {

                                                int iloc = -1;

                                                for (j = 0; j < atom_tmp.size(); ++j) {
                                                    if (atom_tmp[j] == jat) {
                                                        iloc = j;
                                                        break;
                                                    }
                                                }

                                                if (iloc == -1) {
                                                    exit("generate_rotational_constraint", "This cannot happen.");
                                                }

                                                for (j = 0; j < 3; ++j) vec_for_rot[j] = 0.0;

                                                int nsize_equiv = (*iter_cluster).cell.size();

                                                for (j = 0; j < nsize_equiv; ++j) {
                                                    for (int k = 0; k < 3; ++k) {
                                                        vec_for_rot[k] += system->x_image[(*iter_cluster).cell[j][
                                                            iloc]][jat][k];
                                                    }
                                                }

                                                for (j = 0; j < 3; ++j) {
                                                    vec_for_rot[j] /= static_cast<double>(nsize_equiv);
                                                }
                                            }


                                            // mu, nu

                                            interaction_index[order + 1] = 3 * jat + mu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                            sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.
                                                find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if (iter_found != list_found.end()) {
                                                arr_constraint[nparams[order - 1] + (*iter_found).mother]
                                                    += (*iter_found).sign * vec_for_rot[nu];
                                            }

                                            // Exchange mu <--> nu and repeat again.

                                            interaction_index[order + 1] = 3 * jat + nu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                            sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.
                                                find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if (iter_found != list_found.end()) {
                                                arr_constraint[nparams[order - 1] + (*iter_found).mother]
                                                    -= (*iter_found).sign * vec_for_rot[mu];
                                            }
                                        }

                                        for (lambda = 0; lambda < order + 1; ++lambda) {

                                            mu_lambda = interaction_index[lambda] % 3;

                                            for (jcrd = 0; jcrd < 3; ++jcrd) {

                                                for (j = 0; j < order + 1; ++j)
                                                    interaction_tmp[j] = interaction_index[j
                                                    ];

                                                interaction_tmp[lambda] = 3 * interaction_atom[lambda] + jcrd;

                                                levi_factor = 0;

                                                for (j = 0; j < 3; ++j) {
                                                    levi_factor += levi_civita(j, mu, nu) * levi_civita(
                                                        j, mu_lambda, jcrd);
                                                }

                                                if (levi_factor == 0) continue;

                                                sort_tail(order + 1, interaction_tmp);

                                                iter_found = list_found_last.find(FcProperty(order + 1, 1.0,
                                                                                             interaction_tmp, 1));
                                                if (iter_found != list_found_last.end()) {
                                                    arr_constraint[(*iter_found).mother]
                                                        += (*iter_found).sign * static_cast<double>(levi_factor);
                                                }
                                            }
                                        }

                                        if (!is_allzero(nparam_sub, arr_constraint)) {

                                            // A Candidate for another constraint found !
                                            // Add to the appropriate set

                                            if (is_allzero(nparam_sub, arr_constraint, nparams[order - 1])) {
                                                const_rotation_self[order - 1].emplace_back(
                                                    nparams[order - 1], arr_constraint);
                                            } else if (is_allzero(nparams[order - 1], arr_constraint)) {
                                                const_rotation_self[order].emplace_back(
                                                    nparam_sub, arr_constraint, nparams[order - 1]);
                                            } else {
                                                const_rotation_cross[order].emplace_back(nparam_sub, arr_constraint);
                                            }
                                        }

                                    } // nu
                                } // mu

                            } // ixyz

                        } while (g.next());

                    } // direction
                } // icrd
            }

            // Additional constraint for the last order.
            // All IFCs over maxorder-th order are neglected.

            if (order == maxorder - 1 && !exclude_last_R) {

                std::vector<int> interaction_list_now(interaction->interaction_pair[order][i]);
                std::sort(interaction_list_now.begin(), interaction_list_now.end());

                nxyz2 = static_cast<int>(pow(static_cast<double>(3), order + 1));
                allocate(xyzcomponent2, nxyz2, order + 1);
                fcs->get_xyzcomponent(order + 1, xyzcomponent2);

                for (icrd = 0; icrd < 3; ++icrd) {

                    interaction_index[0] = 3 * interaction_atom[0] + icrd;

                    CombinationWithRepetition<int> g_now(interaction_list_now.begin(),
                                                         interaction_list_now.end(), order + 1);
                    do {

                        std::vector<int> data = g_now.now();

                        for (int idata = 0; idata < data.size(); ++idata)
                            interaction_atom[idata + 1] = data[idata];

                        for (ixyz = 0; ixyz < nxyz2; ++ixyz) {

                            for (j = 0; j < order + 1; ++j)
                                interaction_index[j + 1] = 3 * interaction_atom[j + 1] + xyzcomponent2[ixyz][j];

                            for (mu = 0; mu < 3; ++mu) {

                                for (nu = 0; nu < 3; ++nu) {

                                    if (!valid_rotation_axis[mu][nu]) continue;

                                    for (j = 0; j < nparams[order]; ++j) arr_constraint_self[j] = 0.0;

                                    for (lambda = 0; lambda < order + 2; ++lambda) {

                                        mu_lambda = interaction_index[lambda] % 3;

                                        for (jcrd = 0; jcrd < 3; ++jcrd) {

                                            for (j = 0; j < order + 2; ++j)
                                                interaction_tmp[j] = interaction_index[j];

                                            interaction_tmp[lambda] = 3 * interaction_atom[lambda] + jcrd;

                                            levi_factor = 0;
                                            for (j = 0; j < 3; ++j) {
                                                levi_factor += levi_civita(j, mu, nu) * levi_civita(j, mu_lambda, jcrd);
                                            }

                                            if (levi_factor == 0) continue;

                                            sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                                    interaction_tmp, 1));
                                            if (iter_found != list_found.end()) {
                                                arr_constraint_self[(*iter_found).mother]
                                                    += (*iter_found).sign * static_cast<double>(levi_factor);
                                            }
                                        } // jcrd
                                    } // lambda

                                    if (!is_allzero(nparams[order], arr_constraint_self)) {
                                        const_rotation_self[order].emplace_back(
                                            ConstraintClass(nparams[order], arr_constraint_self));
                                    }

                                } // nu
                            } // mu

                        } // ixyz


                    } while (g_now.next());

                } // icrd

                deallocate(xyzcomponent2);
            }
        } // iat

        std::cout << " done." << std::endl;

        if (order > 0) {
            deallocate(xyzcomponent);
        }
        deallocate(arr_constraint);
        deallocate(arr_constraint_self);
        deallocate(interaction_tmp);
        deallocate(interaction_index);
        deallocate(interaction_atom);
    } // order

    for (order = 0; order < maxorder; ++order) {
        nparam_sub = nparams[order] + nparams[order - 1];
        remove_redundant_rows(nparam_sub, const_rotation_cross[order], eps6);
        remove_redundant_rows(nparams[order], const_rotation_self[order], eps6);
    }

    std::cout << "  Finished !" << std::endl << std::endl;

    deallocate(ind);
    deallocate(nparams);
}


int Constraint::levi_civita(const int i, const int j, const int k)
{
    return (j - i) * (k - i) * (k - j) / 2;
}


void Constraint::setup_rotation_axis(bool flag[3][3])
{
    for (auto mu = 0; mu < 3; ++mu) {
        for (auto nu = 0; nu < 3; ++nu) {
            if (mu == nu) {
                flag[mu][nu] = false;
            } else {
                flag[mu][nu] = true;
            }
        }
    }
    std::sort(rotation_axis.begin(), rotation_axis.end());

    if (rotation_axis == "x") {
        flag[0][1] = false;
        flag[1][0] = false;
        flag[0][2] = false;
        flag[2][0] = false;
    } else if (rotation_axis == "y") {
        flag[0][1] = false;
        flag[1][0] = false;
        flag[1][2] = false;
        flag[2][1] = false;
    } else if (rotation_axis == "z") {
        flag[0][2] = false;
        flag[2][0] = false;
        flag[1][2] = false;
        flag[2][1] = false;
    } else if (rotation_axis == "xy") {
        flag[0][1] = false;
        flag[1][0] = false;
    } else if (rotation_axis == "yz") {
        flag[1][2] = false;
        flag[2][1] = false;
    } else if (rotation_axis == "xz") {
        flag[0][2] = false;
        flag[2][0] = false;
    } else if (rotation_axis == "xyz") {
        // do nothing
    } else {
        warn("setup_rotation_axis",
             "Invalid rotation_axis. Default value(xyz) will be used.");
    }
}


void Constraint::fix_forceconstants_to_file(const int order,
                                            Symmetry *symmetry,
                                            Fcs *fcs,
                                            const std::string file_to_fix,
                                            std::vector<ConstraintTypeFix> &const_out)
{
    using namespace boost::property_tree;
    ptree pt;

    try {
        read_xml(file_to_fix, pt);
    }
    catch (std::exception &e) {
        if (order == 0) {
            auto str_error = "Cannot open file FC2XML ( " + file_to_fix + " )";
        } else if (order == 1) {
            auto str_error = "Cannot open file FC3XML ( " + file_to_fix + " )";
        }
        exit("fix_forceconstants_to_file", file_to_fix.c_str());
    }

    int nat_ref = boost::lexical_cast<int>(
        get_value_from_xml(pt, "Data.Structure.NumberOfAtoms"));
    int ntran_ref = boost::lexical_cast<int>(
        get_value_from_xml(pt, "Data.Symmetry.NumberOfTranslations"));
    int natmin_ref = nat_ref / ntran_ref;

    if (natmin_ref != symmetry->nat_prim) {
        exit("fix_forceconstants_to_file",
             "The number of atoms in the primitive cell is not consistent.");
    }

    int nfcs = fcs->nequiv[order].size();

    if (order == 0) {
        int nfcs_ref = boost::lexical_cast<int>(
            get_value_from_xml(pt, "Data.ForceConstants.HarmonicUnique.NFC2"));

        if (nfcs_ref != nfcs) {
            exit("load_reference_system_xml",
                 "The number of harmonic force constants is not consistent.");
        }
    } else if (order == 1) {
        int nfcs_ref = boost::lexical_cast<int>(
            get_value_from_xml(pt, "Data.ForceConstants.CubicUnique.NFC3"));

        if (nfcs_ref != nfcs) {
            exit("load_reference_system_xml",
                 "The number of cubic force constants is not consistent.");
        }
    }

    int **intpair_ref;
    double *fcs_ref;

    allocate(fcs_ref, nfcs);
    allocate(intpair_ref, nfcs, 3);

    int counter = 0;


    if (order == 0) {
        BOOST_FOREACH(const ptree::value_type & child_, pt.get_child("Data.ForceConstants.HarmonicUnique")) {
            if (child_.first == "FC2") {
                const ptree &child = child_.second;
                const std::string str_intpair = child.get<std::string>("<xmlattr>.pairs");
                const std::string str_multiplicity = child.get<std::string>("<xmlattr>.multiplicity");

                std::istringstream is(str_intpair);
                is >> intpair_ref[counter][0] >> intpair_ref[counter][1];
                fcs_ref[counter] = boost::lexical_cast<double>(child.data());
                ++counter;
            }
        }
    } else if (order == 1) {
        BOOST_FOREACH(const ptree::value_type & child_, pt.get_child("Data.ForceConstants.CubicUnique")) {
            if (child_.first == "FC3") {
                const ptree &child = child_.second;
                const std::string str_intpair = child.get<std::string>("<xmlattr>.pairs");
                const std::string str_multiplicity = child.get<std::string>("<xmlattr>.multiplicity");

                std::istringstream is(str_intpair);
                is >> intpair_ref[counter][0] >> intpair_ref[counter][1] >> intpair_ref[counter][2];
                fcs_ref[counter] = boost::lexical_cast<double>(child.data());
                ++counter;
            }
        }
    }

    deallocate(fcs_ref);

    int nterms = order + 2;

    std::unordered_set<FcProperty> list_found;
    std::unordered_set<FcProperty>::iterator iter_found;

    list_found.clear();

    for (auto &list_tmp : fcs->fc_table[order]) {
        list_found.insert(FcProperty(list_tmp));
    }

    for (int i = 0; i < nfcs; ++i) {
        iter_found = list_found.find(FcProperty(nterms, 1.0,
                                                intpair_ref[i], 1));
        if (iter_found == list_found.end()) {
            exit("fix_forceconstants_to_file",
                 "Cannot find equivalent force constant, number: ",
                 i + 1);
        }
        const_out.emplace_back(ConstraintTypeFix((*iter_found).mother, fcs_ref[i]));
    }
    deallocate(intpair_ref);

    list_found.clear();
}


bool Constraint::is_allzero(const int n, const double *arr, const int nshift)
{
    for (int i = nshift; i < n; ++i) {
        if (std::abs(arr[i]) > eps10) {
            return false;
        }
    }
    return true;
}

bool Constraint::is_allzero(const std::vector<int> &vec, int &loc)
{
    loc = -1;
    for (int i = 0; i < vec.size(); ++i) {
        if (std::abs(vec[i]) > 0) {
            loc = i;
            return false;
        }
    }
    return true;
}

bool Constraint::is_allzero(const std::vector<double> &vec, const double tol, int &loc)
{
    loc = -1;
    auto n = vec.size();
    for (int i = 0; i < n; ++i) {
        if (std::abs(vec[i]) > tol) {
            loc = i;
            return false;
        }
    }
    return true;
}


void Constraint::remove_redundant_rows(const int n,
                                       std::vector<ConstraintClass> &Constraint_vec,
                                       const double tolerance)
{
    int i, j;

    int nparam = n;
    int nconst = Constraint_vec.size();
    double *arr_tmp;
    double **mat_tmp;

    int nrank;

    if (nconst > 0) {

        allocate(mat_tmp, nconst, nparam);

        i = 0;

        for (auto &p : Constraint_vec) {
            for (j = 0; j < nparam; ++j) {
                mat_tmp[i][j] = p.w_const[j];
            }
            ++i;
        }

        rref(nconst, nparam, mat_tmp, nrank, tolerance);

        allocate(arr_tmp, nparam);

        Constraint_vec.clear();

        for (i = 0; i < nrank; ++i) {
            for (j = 0; j < i; ++j) arr_tmp[j] = 0.0;

            for (j = i; j < nparam; ++j) {
                if (std::abs(mat_tmp[i][j]) < tolerance) {
                    arr_tmp[j] = 0.0;
                } else {
                    arr_tmp[j] = mat_tmp[i][j];
                }
            }
            Constraint_vec.emplace_back(nparam, arr_tmp);
        }

        deallocate(mat_tmp);
        deallocate(arr_tmp);

    }
}
