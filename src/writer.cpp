/*
 writer.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "writer.h"
#include "alm.h"
#include "constraint.h"
#include "error.h"
#include "fcs.h"
#include "files.h"
#include "fitting.h"
#include "interaction.h"
#include "lasso.h"
#include "memory.h"
#include "patterndisp.h"
#include "symmetry.h"
#include "system.h"
#include "timer.h"
#include "version.h"
#include <iostream>
#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>

using namespace ALM_NS;

Writer::Writer() {}

Writer::~Writer() {}

void Writer::write_input_vars(const ALM *alm) const
{
    unsigned int i;
    int nat, nkd;

    nat = alm->system->get_supercell().number_of_atoms;
    nkd = alm->system->get_supercell().number_of_elems;

    alm->timer->start_clock("writer");

    std::cout << std::endl;
    std::cout << " Input variables:" << std::endl;
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << " General:" << std::endl;
    std::cout << "  PREFIX = " << alm->files->job_title << std::endl;
    std::cout << "  MODE = " << alm->get_run_mode() << std::endl;
    std::cout << "  NAT = " << nat << "; NKD = " << nkd << std::endl;
    std::cout << "  PRINTSYM = " << alm->symmetry->printsymmetry
        << "; TOLERANCE = " << alm->symmetry->tolerance << std::endl;
    std::cout << "  KD = ";
    for (i = 0; i < nkd; ++i) std::cout << std::setw(4) << alm->system->get_kdname()[i];
    std::cout << std::endl;
    std::cout << "  PERIODIC = ";
    for (i = 0; i < 3; ++i) std::cout << std::setw(3) << alm->system->get_periodicity()[i];
    std::cout << std::endl;
    std::cout << "  MAGMOM = " << alm->system->get_str_magmom() << std::endl;
    std::cout << "  HESSIAN = " << alm->files->print_hessian << std::endl;
    std::cout << std::endl;

    std::cout << " Interaction:" << std::endl;
    std::cout << "  NORDER = " << alm->interaction->maxorder << std::endl;
    std::cout << "  NBODY = ";
    for (i = 0; i < alm->interaction->maxorder; ++i)
        std::cout << std::setw(3) << alm->interaction->nbody_include[i];

    std::cout << std::endl << std::endl;


    if (alm->get_run_mode() == "suggest") {
        std::cout << "  DBASIS = " << alm->displace->disp_basis << std::endl;
        std::cout << std::endl;

    } else if (alm->get_run_mode() == "fitting") {
        std::cout << " Fitting:" << std::endl;
        std::cout << "  DFILE = " << alm->files->file_disp << std::endl;
        std::cout << "  FFILE = " << alm->files->file_force << std::endl;
        std::cout << "  NDATA = " << alm->fitting->ndata << "; NSTART = " << alm->fitting->nstart
            << "; NEND = " << alm->fitting->nend << std::endl;
        std::cout << "  ICONST = " << alm->constraint->constraint_mode << std::endl;
        std::cout << "  ROTAXIS = " << alm->constraint->rotation_axis << std::endl;
        std::cout << "  FC2XML = " << alm->constraint->fc2_file << std::endl;
        std::cout << "  FC3XML = " << alm->constraint->fc3_file << std::endl;
        std::cout << "  SPARSE = " << alm->fitting->use_sparseQR << std::endl;
        std::cout << std::endl;
    } else if (alm->get_run_mode() == "lasso") {
        std::cout << " Fitting:" << std::endl;
        std::cout << "  DFILE = " << alm->files->file_disp << std::endl;
        std::cout << "  FFILE = " << alm->files->file_force << std::endl;
        std::cout << "  NDATA = " << alm->fitting->ndata << "; NSTART = " << alm->fitting->nstart
            << "; NEND = " << alm->fitting->nend << std::endl;
        std::cout << "  SKIP = " << alm->fitting->skip_s + 1 << "-" << alm->fitting->skip_e << std::endl;
        std::cout << "  ICONST = " << alm->constraint->constraint_mode << std::endl;
        std::cout << "  ROTAXIS = " << alm->constraint->rotation_axis << std::endl;
        std::cout << "  FC2XML = " << alm->constraint->fc2_file << std::endl;
        std::cout << "  FC3XML = " << alm->constraint->fc3_file << std::endl;
        std::cout << std::endl;
        std::cout << " Lasso:" << std::endl;
        std::cout << "  LASSO_ALPHA = " << alm->lasso->l1_alpha << std::endl;
        std::cout << "  LASSO_MINALPHA = " << alm->lasso->l1_alpha_min;
        std::cout << "; LASSO_MAXALPHA = " << alm->lasso->l1_alpha_max << std::endl;
        std::cout << "  LASSO_NALPHA = " << alm->lasso->num_l1_alpha << std::endl;
        std::cout << "  STANDARDIZE = " << alm->lasso->standardize << std::endl;
        std::cout << "  LASSO_DNORM = " << alm->lasso->disp_norm << std::endl;
        std::cout << "  LASSO_TOL = " << alm->lasso->lasso_tol << std::endl;
        std::cout << "  LASSO_MAXITER = " << alm->lasso->maxiter << std::endl;
        std::cout << "  LASSO_CV = " << std::setw(5) << alm->lasso->lasso_cv << std::endl;
        std::cout << "  LASSO_FREQ = " << std::setw(5) << alm->lasso->output_frequency << std::endl;
        std::cout << "  LASSO_ALGO = " << std::setw(5) << alm->lasso->lasso_algo << std::endl;
        if (alm->lasso->lasso_algo == 1) {
            std::cout << "  LASSO_LAMBDA = " << alm->lasso->l2_lambda << std::endl;
            std::cout << "  LASSO_ZERO_THR = " << std::setw(15) << alm->lasso->lasso_zero_thr << std::endl;
            std::cout << "  LASSO_MAXITER_CG = " << std::setw(5) << alm->lasso->maxiter_cg << std::endl;
            std::cout << "  LASSO_PCG = " << std::setw(5) << alm->lasso->lasso_pcg << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    alm->timer->stop_clock("writer");
}

void Writer::writeall(ALM *alm)
{
    alm->timer->start_clock("writer");

    if (alm->get_verbosity() > 0)
        std::cout << " The following files are created:" << std::endl << std::endl;

    write_force_constants(alm);
    // write_misc_xml breaks data in fcs.
    write_misc_xml(alm);
    if (alm->files->print_hessian) write_hessian(alm);
    //   write_in_QEformat(alm);

    alm->timer->stop_clock("writer");
}

void Writer::write_force_constants(ALM *alm) const
{
    int order, j, l, m;
    unsigned int ui;
    int multiplicity;
    double distmax;
    std::string *str_fcs;
    std::string str_tmp;
    std::ofstream ofs_fcs;
    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<InteractionCluster>::iterator iter_cluster;

    int maxorder = alm->interaction->maxorder;

    ofs_fcs.open(alm->files->file_fcs.c_str(), std::ios::out);
    if (!ofs_fcs) exit("write_force_constants", "cannot open fcs file");

    ofs_fcs << " *********************** Force Constants (FCs) ***********************" << std::endl;
    ofs_fcs << " *        Force constants are printed in Rydberg atomic units.       *" << std::endl;
    ofs_fcs << " *        FC2: Ry/a0^2     FC3: Ry/a0^3     FC4: Ry/a0^4   etc.      *" << std::endl;
    ofs_fcs << " *        FC?: Ry/a0^?     a0 = Bohr radius                          *" << std::endl;
    ofs_fcs << " *                                                                   *" << std::endl;
    ofs_fcs << " *        The value shown in the last column is the distance         *" << std::endl;
    ofs_fcs << " *        between the most distant atomic pairs.                     *" << std::endl;
    ofs_fcs << " *********************************************************************" << std::endl;
    ofs_fcs << std::endl;
    ofs_fcs << " ----------------------------------------------------------------------" << std::endl;
    ofs_fcs << "      Index              FCs         P        Pairs     Distance [Bohr]" << std::endl;
    ofs_fcs << " (Global, Local)              (Multiplicity)                           " << std::endl;
    ofs_fcs << " ----------------------------------------------------------------------" << std::endl;

    allocate(str_fcs, maxorder);

    for (order = 0; order < maxorder; ++order) {
        str_fcs[order] = "*FC" + std::to_string(order + 2);
    }

    int k = 0;

    for (order = 0; order < maxorder; ++order) {

        m = 0;

        if (!alm->fcs->nequiv[order].empty()) {

            ofs_fcs << std::endl << std::setw(6) << str_fcs[order] << std::endl;

            for (ui = 0; ui < alm->fcs->nequiv[order].size(); ++ui) {

                ofs_fcs << std::setw(8) << k + 1 << std::setw(8) << ui + 1
                    << std::setw(18) << std::setprecision(7)
                    << std::scientific << alm->fitting->params[k];

                atom_tmp.clear();
                for (l = 1; l < order + 2; ++l) {
                    atom_tmp.push_back(alm->fcs->fc_table[order][m].elems[l] / 3);
                }
                j = alm->symmetry->map_s2p[alm->fcs->fc_table[order][m].elems[0] / 3].atom_num;
                std::sort(atom_tmp.begin(), atom_tmp.end());

                iter_cluster = alm->interaction->interaction_cluster[order][j].find(
                    InteractionCluster(atom_tmp, cell_dummy));

                if (iter_cluster == alm->interaction->interaction_cluster[order][j].end()) {
                    std::cout << std::setw(5) << j;
                    for (l = 0; l < order + 1; ++l) {
                        std::cout << std::setw(5) << atom_tmp[l];
                    }
                    std::cout << std::endl;
                    exit("write_force_constants",
                         "This cannot happen.");
                }

                multiplicity = (*iter_cluster).cell.size();
                distmax = (*iter_cluster).distmax;
                ofs_fcs << std::setw(4) << multiplicity;

                for (l = 0; l < order + 2; ++l) {
                    ofs_fcs << std::setw(7)
                        << easyvizint(alm->fcs->fc_table[order][m].elems[l]);
                }
                ofs_fcs << std::setw(12) << std::setprecision(3)
                    << std::fixed << distmax << std::endl;

                m += alm->fcs->nequiv[order][ui];
                ++k;
            }
        }
    }

    ofs_fcs << std::endl;

    /*  if (alm->constraint->extra_constraint_from_symmetry) {

          ofs_fcs << " -------------- Constraints from crystal symmetry --------------" << std::endl << std::endl;
          for (order = 0; order < maxorder; ++order) {
              int nparam = alm->fcs->nequiv[order].size();


              for (auto p = alm->constraint->const_symmetry[order].begin();
                   p != alm->constraint->const_symmetry[order].end();
                   ++p) {
                  ofs_fcs << "   0 = " << std::scientific << std::setprecision(6);
                  auto const_pointer = *p;
                  for (j = 0; j < nparam; ++j) {
                      if (std::abs(const_pointer.w_const[j]) > eps8) {
                          str_tmp = " * (FC" + std::to_string(order + 2)
                              + "_" + std::to_string(j + 1) + ")";
                          ofs_fcs << std::setw(10) << std::right
                              << std::showpos << const_pointer.w_const[j];
                          ofs_fcs << std::setw(12) << std::left << str_tmp;
                      }
                  }
                  ofs_fcs << std::endl;
              }
              ofs_fcs << std::endl;
          }
          ofs_fcs << std::endl;
      }*/

    ofs_fcs.unsetf(std::ios::showpos);

    for (order = 0; order < maxorder; ++order) {
        str_fcs[order] = "**FC" + std::to_string(order + 2);
    }

    ofs_fcs << std::endl << std::endl;
    ofs_fcs << " ------------------------ All FCs below ------------------------" << std::endl;

    int ip = 0;
    int id;

    for (order = 0; order < maxorder; ++order) {

        id = 0;

        if (!alm->fcs->nequiv[order].empty()) {
            ofs_fcs << std::endl << std::setw(6) << str_fcs[order] << std::endl;

            for (unsigned int iuniq = 0; iuniq < alm->fcs->nequiv[order].size(); ++iuniq) {

                str_tmp = "  # FC" + std::to_string(order + 2) + "_";
                str_tmp += std::to_string(iuniq + 1);

                ofs_fcs << str_tmp << std::setw(5) << alm->fcs->nequiv[order][iuniq]
                    << std::setw(16) << std::scientific
                    << std::setprecision(7) << alm->fitting->params[ip] << std::endl;

                for (j = 0; j < alm->fcs->nequiv[order][iuniq]; ++j) {
                    ofs_fcs << std::setw(5) << j + 1 << std::setw(12)
                        << std::setprecision(5) << std::fixed << alm->fcs->fc_table[order][id].sign;
                    for (k = 0; k < order + 2; ++k) {
                        ofs_fcs << std::setw(6)
                            << easyvizint(alm->fcs->fc_table[order][id].elems[k]);
                    }
                    ofs_fcs << std::endl;
                    ++id;
                }
                ofs_fcs << std::endl;
                ++ip;
            }
        }
    }
    deallocate(str_fcs);
    ofs_fcs.close();

    if (alm->get_verbosity() > 0) {
        std::cout << " Force constants in a human-readable format : "
            << alm->files->file_fcs << std::endl;
    }
}

void Writer::write_displacement_pattern(ALM *alm) const
{
    int maxorder = alm->interaction->maxorder;

    std::ofstream ofs_pattern;
    std::string file_disp_pattern;

    if (alm->get_verbosity() > 0)
        std::cout << " Suggested displacement patterns are printed in the following files: " << std::endl;

    for (int order = 0; order < maxorder; ++order) {

        if (order == 0) {
            file_disp_pattern = alm->files->job_title + ".pattern_HARMONIC";
        } else {
            file_disp_pattern = alm->files->job_title + ".pattern_ANHARM"
                + std::to_string(order + 2);
        }

        ofs_pattern.open(file_disp_pattern.c_str(), std::ios::out);
        if (!ofs_pattern) {
            exit("write_displacement_pattern",
                 "Cannot open file_disp_pattern");
        }

        int counter = 0;

        ofs_pattern << "Basis : " << alm->displace->disp_basis[0] << std::endl;

        for (auto entry : alm->displace->pattern_all[order]) {
            ++counter;

            ofs_pattern << std::setw(5) << counter << ":"
                << std::setw(5) << entry.atoms.size() << std::endl;
            for (int i = 0; i < entry.atoms.size(); ++i) {
                ofs_pattern << std::setw(7) << entry.atoms[i] + 1;
                for (int j = 0; j < 3; ++j) {
                    ofs_pattern << std::setw(15) << entry.directions[3 * i + j];
                }
                ofs_pattern << std::endl;
            }
        }

        ofs_pattern.close();

        if (alm->get_verbosity() > 0) {
            std::cout << "  " << alm->interaction->str_order[order]
                << " : " << file_disp_pattern << std::endl;
        }

    }
    if (alm->get_verbosity() > 0) std::cout << std::endl;
}


void Writer::write_misc_xml(ALM *alm)
{
    SystemInfo system_structure;

    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            system_structure.lattice_vector[i][j]
                = alm->system->get_supercell().lattice_vector[i][j];
        }
    }

    system_structure.nat = alm->system->get_supercell().number_of_atoms;
    system_structure.natmin = alm->symmetry->nat_prim;
    system_structure.ntran = alm->symmetry->ntran;
    system_structure.nspecies = alm->system->get_supercell().number_of_elems;

    AtomProperty prop_tmp;

    for (i = 0; i < alm->system->get_supercell().number_of_atoms; ++i) {
        prop_tmp.x = alm->system->get_supercell().x_fractional[i][0];
        prop_tmp.y = alm->system->get_supercell().x_fractional[i][1];
        prop_tmp.z = alm->system->get_supercell().x_fractional[i][2];
        prop_tmp.kind = alm->system->get_supercell().kind[i];
        prop_tmp.atom = alm->symmetry->map_s2p[i].atom_num + 1;
        prop_tmp.tran = alm->symmetry->map_s2p[i].tran_num + 1;

        system_structure.atoms.emplace_back(AtomProperty(prop_tmp));
    }

    using boost::property_tree::ptree;

    ptree pt;
    std::string str_pos[3];

    pt.put("Data.ALM_version", ALAMODE_VERSION);
    pt.put("Data.Fitting.DisplaceFile", alm->files->file_disp);
    pt.put("Data.Fitting.ForceFile", alm->files->file_force);
    pt.put("Data.Fitting.Constraint", alm->constraint->constraint_mode);

    pt.put("Data.Structure.NumberOfAtoms", system_structure.nat);
    pt.put("Data.Structure.NumberOfElements", system_structure.nspecies);

    for (i = 0; i < system_structure.nspecies; ++i) {
        ptree &child = pt.add("Data.Structure.AtomicElements.element",
                              alm->system->get_kdname()[i]);
        child.put("<xmlattr>.number", i + 1);
    }

    for (i = 0; i < 3; ++i) {
        str_pos[i].clear();
        for (j = 0; j < 3; ++j) {
            str_pos[i] += " " + double2string(system_structure.lattice_vector[j][i]);
        }
    }
    pt.put("Data.Structure.LatticeVector", "");
    pt.put("Data.Structure.LatticeVector.a1", str_pos[0]);
    pt.put("Data.Structure.LatticeVector.a2", str_pos[1]);
    pt.put("Data.Structure.LatticeVector.a3", str_pos[2]);

    std::stringstream ss;
    ss << alm->system->get_periodicity()[0] << " "
        << alm->system->get_periodicity()[1] << " "
        << alm->system->get_periodicity()[2];
    pt.put("Data.Structure.Periodicity", ss.str());

    pt.put("Data.Structure.Position", "");
    std::string str_tmp;

    for (i = 0; i < system_structure.nat; ++i) {
        str_tmp.clear();
        for (j = 0; j < 3; ++j) str_tmp += " " + double2string(alm->system->get_supercell().x_fractional[i][j]);
        ptree &child = pt.add("Data.Structure.Position.pos", str_tmp);
        child.put("<xmlattr>.index", i + 1);
        child.put("<xmlattr>.element", alm->system->get_kdname()[alm->system->get_supercell().kind[i] - 1]);
    }

    pt.put("Data.Symmetry.NumberOfTranslations", alm->symmetry->ntran);
    for (i = 0; i < system_structure.ntran; ++i) {
        for (j = 0; j < system_structure.natmin; ++j) {
            ptree &child = pt.add("Data.Symmetry.Translations.map",
                                  alm->symmetry->map_p2s[j][i] + 1);
            child.put("<xmlattr>.tran", i + 1);
            child.put("<xmlattr>.atom", j + 1);
        }
    }

    if (alm->system->get_spin().lspin) {
        pt.put("Data.MagneticMoments", "");
        pt.put("Data.MagneticMoments.Noncollinear", alm->system->get_spin().noncollinear);
        pt.put("Data.MagneticMoments.TimeReversalSymmetry", alm->system->get_spin().time_reversal_symm);
        for (i = 0; i < system_structure.nat; ++i) {
            str_tmp.clear();
            for (j = 0; j < 3; ++j) str_tmp += " " + double2string(alm->system->get_spin().magmom[i][j], 5);
            ptree &child = pt.add("Data.MagneticMoments.mag", str_tmp);
            child.put("<xmlattr>.index", i + 1);
        }
    }

    pt.put("Data.ForceConstants", "");
    str_tmp.clear();

    pt.put("Data.ForceConstants.HarmonicUnique.NFC2", alm->fcs->nequiv[0].size());

    int ihead = 0;
    int k = 0;
    int nelem = alm->interaction->maxorder + 1;
    int *pair_tmp;
    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<InteractionCluster>::iterator iter_cluster;
    int multiplicity;


    allocate(pair_tmp, nelem);

    for (unsigned int ui = 0; ui < alm->fcs->nequiv[0].size(); ++ui) {

        for (i = 0; i < 2; ++i) {
            pair_tmp[i] = alm->fcs->fc_table[0][ihead].elems[i] / 3;
        }
        j = alm->symmetry->map_s2p[pair_tmp[0]].atom_num;

        atom_tmp.clear();
        atom_tmp.push_back(pair_tmp[1]);

        iter_cluster = alm->interaction->interaction_cluster[0][j].find(
            InteractionCluster(atom_tmp, cell_dummy));
        if (iter_cluster == alm->interaction->interaction_cluster[0][j].end()) {
            exit("load_reference_system_xml",
                 "Cubic force constant is not found.");
        }

        multiplicity = (*iter_cluster).cell.size();

        ptree &child = pt.add("Data.ForceConstants.HarmonicUnique.FC2",
                              double2string(alm->fitting->params[k]));
        child.put("<xmlattr>.pairs",
                  std::to_string(alm->fcs->fc_table[0][ihead].elems[0])
                  + " " + std::to_string(alm->fcs->fc_table[0][ihead].elems[1]));
        child.put("<xmlattr>.multiplicity", multiplicity);
        ihead += alm->fcs->nequiv[0][ui];
        ++k;
    }
    ihead = 0;


    if (alm->interaction->maxorder > 1) {

        pt.put("Data.ForceConstants.CubicUnique.NFC3", alm->fcs->nequiv[1].size());

        for (unsigned int ui = 0; ui < alm->fcs->nequiv[1].size(); ++ui) {
            for (i = 0; i < 3; ++i) {
                pair_tmp[i] = alm->fcs->fc_table[1][ihead].elems[i] / 3;
            }
            j = alm->symmetry->map_s2p[pair_tmp[0]].atom_num;

            atom_tmp.clear();
            for (i = 1; i < 3; ++i) {
                atom_tmp.push_back(pair_tmp[i]);
            }
            std::sort(atom_tmp.begin(), atom_tmp.end());

            iter_cluster = alm->interaction->interaction_cluster[1][j].find(
                InteractionCluster(atom_tmp, cell_dummy));
            if (iter_cluster == alm->interaction->interaction_cluster[1][j].end()) {
                exit("load_reference_system_xml",
                     "Cubic force constant is not found.");
            }
            multiplicity = (*iter_cluster).cell.size();


            ptree &child = pt.add("Data.ForceConstants.CubicUnique.FC3",
                                  double2string(alm->fitting->params[k]));
            child.put("<xmlattr>.pairs",
                      std::to_string(alm->fcs->fc_table[1][ihead].elems[0])
                      + " " + std::to_string(alm->fcs->fc_table[1][ihead].elems[1])
                      + " " + std::to_string(alm->fcs->fc_table[1][ihead].elems[2]));
            child.put("<xmlattr>.multiplicity", multiplicity);
            ihead += alm->fcs->nequiv[1][ui];
            ++k;
        }
    }

    int ip;
    int imult;
    std::string elementname = "Data.ForceConstants.HARMONIC.FC2";

    std::sort(alm->fcs->fc_table[0].begin(), alm->fcs->fc_table[0].end());

    for (auto it = alm->fcs->fc_table[0].begin(); it != alm->fcs->fc_table[0].end(); ++it) {
        FcProperty fctmp = *it;
        ip = fctmp.mother;

        for (k = 0; k < 2; ++k) {
            pair_tmp[k] = fctmp.elems[k] / 3;
        }

        j = alm->symmetry->map_s2p[pair_tmp[0]].atom_num;

        atom_tmp.clear();
        atom_tmp.push_back(pair_tmp[1]);

        iter_cluster = alm->interaction->interaction_cluster[0][j].find(
            InteractionCluster(atom_tmp, cell_dummy));

        if (iter_cluster != alm->interaction->interaction_cluster[0][j].end()) {
            multiplicity = (*iter_cluster).cell.size();

            for (imult = 0; imult < multiplicity; ++imult) {
                std::vector<int> cell_now = (*iter_cluster).cell[imult];

                ptree &child = pt.add(elementname,
                                      double2string(alm->fitting->params[ip] * fctmp.sign
                                          / static_cast<double>(multiplicity)));

                child.put("<xmlattr>.pair1", std::to_string(j + 1)
                          + " " + std::to_string(fctmp.elems[0] % 3 + 1));

                for (k = 1; k < 2; ++k) {
                    child.put("<xmlattr>.pair" + std::to_string(k + 1),
                              std::to_string(pair_tmp[k] + 1)
                              + " " + std::to_string(fctmp.elems[k] % 3 + 1)
                              + " " + std::to_string(cell_now[k - 1] + 1));
                }
            }
        } else {
            exit("write_misc_xml", "This cannot happen.");
        }
    }

    int ishift = alm->fcs->nequiv[0].size();

    // Print anharmonic force constants to the xml file.


    int order;
    for (order = 1; order < alm->interaction->maxorder; ++order) {

        std::sort(alm->fcs->fc_table[order].begin(), alm->fcs->fc_table[order].end());

        for (auto it = alm->fcs->fc_table[order].begin();
             it != alm->fcs->fc_table[order].end(); ++it) {
            FcProperty fctmp = *it;
            ip = fctmp.mother + ishift;

            for (k = 0; k < order + 2; ++k) {
                pair_tmp[k] = fctmp.elems[k] / 3;
            }
            j = alm->symmetry->map_s2p[pair_tmp[0]].atom_num;

            atom_tmp.clear();

            for (k = 1; k < order + 2; ++k) {
                atom_tmp.push_back(pair_tmp[k]);
            }
            std::sort(atom_tmp.begin(), atom_tmp.end());

            elementname = "Data.ForceConstants.ANHARM"
                + std::to_string(order + 2)
                + ".FC" + std::to_string(order + 2);


            iter_cluster = alm->interaction->interaction_cluster[order][j].find(
                InteractionCluster(atom_tmp, cell_dummy));

            if (iter_cluster != alm->interaction->interaction_cluster[order][j].end()) {
                multiplicity = (*iter_cluster).cell.size();

                for (imult = 0; imult < multiplicity; ++imult) {
                    std::vector<int> cell_now = (*iter_cluster).cell[imult];

                    ptree &child = pt.add(elementname,
                                          double2string(alm->fitting->params[ip] * fctmp.sign
                                              / static_cast<double>(multiplicity)));

                    child.put("<xmlattr>.pair1", std::to_string(j + 1)
                              + " " + std::to_string(fctmp.elems[0] % 3 + 1));

                    for (k = 1; k < order + 2; ++k) {
                        child.put("<xmlattr>.pair" + std::to_string(k + 1),
                                  std::to_string(pair_tmp[k] + 1)
                                  + " " + std::to_string(fctmp.elems[k] % 3 + 1)
                                  + " " + std::to_string(cell_now[k - 1] + 1));
                    }
                }
            } else {
                exit("write_misc_xml", "This cannot happen.");
            }
        }
        ishift += alm->fcs->nequiv[order].size();
    }

    using namespace boost::property_tree::xml_parser;
    const int indent = 2;

    std::string file_xml = alm->files->job_title + ".xml";

#if BOOST_VERSION >= 105600
    write_xml(file_xml, pt, std::locale(),
              xml_writer_make_settings<ptree::key_type>(' ', indent,
                                                        widen<std::string>("utf-8")));
#else
    write_xml(file_xml, pt, std::locale(),
              xml_writer_make_settings(' ', indent, widen<char>("utf-8")));
#endif

    deallocate(pair_tmp);

    if (alm->get_verbosity() > 0) {
        std::cout << " Input data for the phonon code ANPHON      : " << file_xml << std::endl;
    }
}

void Writer::write_hessian(ALM *alm) const
{
    int i, j, itran, ip;
    int pair_tmp[2];
    int pair_tran[2];
    std::ofstream ofs_hes;
    double **hessian;

    //ALMCore *alm = alm->get_alm();
    int nat3 = 3 * alm->system->get_supercell().number_of_atoms;

    allocate(hessian, nat3, nat3);

    for (i = 0; i < nat3; ++i) {
        for (j = 0; j < nat3; ++j) {
            hessian[i][j] = 0.0;
        }
    }

    for (auto it = alm->fcs->fc_table[0].begin();
         it != alm->fcs->fc_table[0].end(); ++it) {
        FcProperty fctmp = *it;
        ip = fctmp.mother;

        for (i = 0; i < 2; ++i) pair_tmp[i] = fctmp.elems[i] / 3;
        for (itran = 0; itran < alm->symmetry->ntran; ++itran) {
            for (i = 0; i < 2; ++i) {
                pair_tran[i] = alm->symmetry->map_sym[pair_tmp[i]][alm->symmetry->symnum_tran[itran]];
            }
            hessian[3 * pair_tran[0] + fctmp.elems[0] % 3][3 * pair_tran[1] + fctmp.elems[1] % 3]
                = alm->fitting->params[ip] * fctmp.sign;
        }
    }

    ofs_hes.open(alm->files->file_hes.c_str(), std::ios::out);
    if (!ofs_hes) exit("write_hessian", "cannot create hessian file");

    ofs_hes << "# atom1, xyz1, atom2, xyz2, FC2 (Ryd/Bohr^2)" << std::endl;
    for (i = 0; i < nat3; ++i) {
        for (j = 0; j < nat3; ++j) {
            ofs_hes << std::setw(5) << i / 3 + 1;
            ofs_hes << std::setw(5) << i % 3 + 1;
            ofs_hes << std::setw(6) << j / 3 + 1;
            ofs_hes << std::setw(5) << j % 3 + 1;
            ofs_hes << std::setw(25) << std::setprecision(15)
                << std::scientific << hessian[i][j];
            ofs_hes << std::endl;
        }
    }
    ofs_hes.close();
    deallocate(hessian);

    if (alm->get_verbosity()) {
        std::cout << " Complete Hessian matrix                    : " << alm->files->file_hes << std::endl;
    }
}

std::string Writer::double2string(const double d,
                                  const int nprec) const
{
    std::string rt;
    std::stringstream ss;

    ss << std::scientific << std::setprecision(nprec) << d;
    ss >> rt;
    return rt;
}

void Writer::write_in_QEformat(ALM *alm) const
{
    int i, j, itran, ip;
    int pair_tmp[2];
    int pair_tran[2];
    std::ofstream ofs_hes;
    double **hessian;
    int nat3 = 3 * alm->system->get_supercell().number_of_atoms;

    allocate(hessian, nat3, nat3);

    for (i = 0; i < nat3; ++i) {
        for (j = 0; j < nat3; ++j) {
            hessian[i][j] = 0.0;
        }
    }

    for (auto it = alm->fcs->fc_table[0].begin(); it != alm->fcs->fc_table[0].end(); ++it) {
        FcProperty fctmp = *it;
        ip = fctmp.mother;

        for (i = 0; i < 2; ++i) pair_tmp[i] = fctmp.elems[i] / 3;
        for (itran = 0; itran < alm->symmetry->ntran; ++itran) {
            for (i = 0; i < 2; ++i) {
                pair_tran[i] = alm->symmetry->map_sym[pair_tmp[i]][alm->symmetry->symnum_tran[itran]];
            }
            hessian[3 * pair_tran[0] + fctmp.elems[0] % 3][3 * pair_tran[1] + fctmp.elems[1] % 3]
                = alm->fitting->params[ip] * fctmp.sign;
        }
    }

    std::string file_fc = alm->files->job_title + ".fc";

    ofs_hes.open(file_fc.c_str(), std::ios::out);
    if (!ofs_hes) exit("write_hessian", "cannot create hessian file");

    ofs_hes << "  1  1  1" << std::endl;
    for (int icrd = 0; icrd < 3; ++icrd) {
        for (int jcrd = 0; jcrd < 3; ++jcrd) {
            for (i = 0; i < alm->system->get_supercell().number_of_atoms; ++i) {
                for (j = 0; j < alm->system->get_supercell().number_of_atoms; ++j) {
                    ofs_hes << std::setw(3) << icrd + 1;
                    ofs_hes << std::setw(3) << jcrd + 1;
                    ofs_hes << std::setw(3) << i + 1;
                    ofs_hes << std::setw(3) << j + 1;
                    ofs_hes << std::endl;
                    ofs_hes << "  1  1  1 " << std::setw(20) << std::setprecision(13)
                        << std::scientific << hessian[3 * j + jcrd][3 * i + icrd];
                    ofs_hes << std::endl;
                }
            }
        }
    }
    ofs_hes.close();
    deallocate(hessian);
}

std::string Writer::easyvizint(const int n) const
{
    int atmn = n / 3 + 1;
    int crdn = n % 3;
    std::string str_crd[3] = {"x", "y", "z"};
    std::string str_tmp = std::to_string(atmn);
    str_tmp += str_crd[crdn];

    return str_tmp;
}
