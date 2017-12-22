/*
 writer.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include "writer.h"
#include "alm.h"
#include "alm_core.h"
#include "system.h"
#include "interaction.h"
#include "memory.h"
#include "symmetry.h"
#include "error.h"
#include "files.h"
#include "fcs.h"
#include "fitting.h"
#include "constraint.h"
#include "patterndisp.h"
#include "version.h"
#include "timer.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>

using namespace ALM_NS;

Writer::Writer()
{
}

Writer::~Writer()
{
}

void Writer::write_input_vars(ALM *alm)
{
    unsigned int i;

    ALMCore *alm_core = alm->get_alm_core();

    alm_core->timer->start_clock("writer");

    std::cout << std::endl;
    std::cout << " Input variables:" << std::endl;
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << " General:" << std::endl;
    std::cout << "  PREFIX = " << alm_core->files->job_title << std::endl;
    std::cout << "  MODE = " << alm_core->mode << std::endl;
    std::cout << "  NAT = " << alm_core->system->nat << "; NKD = " << alm_core->system->nkd << std::endl;
    std::cout << "  NSYM = " << alm_core->symmetry->nsym << "; PRINTSYM = " << alm_core->symmetry->printsymmetry
        << "; TOLERANCE = " << alm_core->symmetry->tolerance << std::endl;
    std::cout << "  KD = ";
    for (i = 0; i < alm_core->system->nkd; ++i) std::cout << std::setw(4) << alm_core->system->kdname[i];
    std::cout << std::endl;
    std::cout << "  PERIODIC = ";
    for (i = 0; i < 3; ++i) std::cout << std::setw(3) << alm_core->interaction->is_periodic[i];
    std::cout << std::endl;
    std::cout << "  MAGMOM = " << alm_core->system->str_magmom << std::endl;
    std::cout << "  HESSIAN = " << alm_core->files->print_hessian << std::endl;
    std::cout << std::endl;


    std::cout << " Interaction:" << std::endl;
    std::cout << "  NORDER = " << alm_core->interaction->maxorder << std::endl;
    std::cout << "  NBODY = ";
    for (i = 0; i < alm_core->interaction->maxorder; ++i)
        std::cout << std::setw(3) << alm_core->interaction->nbody_include[i];

    std::cout << std::endl << std::endl;


    if (alm_core->mode == "suggest") {
        std::cout << "  DBASIS = " << alm_core->displace->disp_basis << std::endl;
        std::cout << std::endl;

    } else if (alm_core->mode == "fitting") {
        std::cout << " Fitting:" << std::endl;
        std::cout << "  DFILE = " << alm_core->files->file_disp << std::endl;
        std::cout << "  FFILE = " << alm_core->files->file_force << std::endl;
        std::cout << "  NDATA = " << alm_core->system->ndata << "; NSTART = " << alm_core->system->nstart
            << "; NEND = " << alm_core->system->nend << std::endl;
        std::cout << "  MULTDAT = " << alm_core->symmetry->multiply_data << std::endl;
        std::cout << "  ICONST = " << alm_core->constraint->constraint_mode << std::endl;
        std::cout << "  ROTAXIS = " << alm_core->constraint->rotation_axis << std::endl;
        std::cout << "  FC2XML = " << alm_core->constraint->fc2_file << std::endl;
        std::cout << "  FC3XML = " << alm_core->constraint->fc3_file << std::endl;
        std::cout << std::endl;
    }
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    alm_core->timer->stop_clock("writer");
}

void Writer::writeall(ALM *alm)
{
    ALMCore *alm_core = alm->get_alm_core();

 //   alm_core->timer->start_clock("writer");

    std::cout << " The following files are created:" << std::endl << std::endl;
    write_force_constants(alm);
    // write_misc_xml breaks data in fcs.
    write_misc_xml(alm);
    if (alm_core->files->print_hessian) write_hessian(alm);
 //   write_in_QEformat(alm);
    std::cout << std::endl;

 //   alm_core->timer->stop_clock("writer");
}

void Writer::write_force_constants(ALM *alm)
{
    int order, j, k, l, m;
    unsigned int ui;
    int multiplicity;
    double distmax;
    std::string *str_fcs;
    std::string str_tmp;
    std::ofstream ofs_fcs;
    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<MinimumDistanceCluster>::iterator iter_cluster;

    ALMCore *alm_core = alm->get_alm_core();
    int maxorder = alm_core->interaction->maxorder;

    ofs_fcs.open(alm_core->files->file_fcs.c_str(), std::ios::out);
    if (!ofs_fcs) alm_core->error->exit("openfiles", "cannot open fcs file");

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
        str_fcs[order] = "*FC" + boost::lexical_cast<std::string>(order + 2);
    }

    k = 0;

    for (order = 0; order < maxorder; ++order) {

        m = 0;

        if (alm_core->fcs->ndup[order].size() > 0) {

            ofs_fcs << std::endl << std::setw(6) << str_fcs[order] << std::endl;

            for (ui = 0; ui < alm_core->fcs->ndup[order].size(); ++ui) {

                ofs_fcs << std::setw(8) << k + 1 << std::setw(8) << ui + 1
                    << std::setw(18) << std::setprecision(7)
                    << std::scientific << alm_core->fitting->params[k];

                atom_tmp.clear();
                for (l = 1; l < order + 2; ++l) {
                    atom_tmp.push_back(alm_core->fcs->fc_table[order][m].elems[l] / 3);
                }
                j = alm_core->symmetry->map_s2p[alm_core->fcs->fc_table[order][m].elems[0] / 3].atom_num;
                std::sort(atom_tmp.begin(), atom_tmp.end());

                iter_cluster = alm_core->interaction->mindist_cluster[order][j].find(
                    MinimumDistanceCluster(atom_tmp, cell_dummy));

                if (iter_cluster != alm_core->interaction->mindist_cluster[order][j].end()) {
                    multiplicity = (*iter_cluster).cell.size();
                    distmax = (*iter_cluster).distmax;
                } else {
                    std::cout << std::setw(5) << j;
                    for (l = 0; l < order + 1; ++l) {
                        std::cout << std::setw(5) << atom_tmp[l];
                    }
                    std::cout << std::endl;
                    alm_core->error->exit("write_force_constants",
                                          "This cannot happen.");
                }
                ofs_fcs << std::setw(4) << multiplicity;

                for (l = 0; l < order + 2; ++l) {
                    ofs_fcs << std::setw(7)
                        << alm_core->fcs->easyvizint(alm_core->fcs->fc_table[order][m].elems[l]);
                }
                ofs_fcs << std::setw(12) << std::setprecision(3)
                    << std::fixed << distmax << std::endl;

                m += alm_core->fcs->ndup[order][ui];
                ++k;
            }
        }
    }

    ofs_fcs << std::endl;

    if (alm_core->constraint->extra_constraint_from_symmetry) {

        ofs_fcs << " -------------- Constraints from crystal symmetry --------------" << std::endl << std::endl;;
        for (order = 0; order < maxorder; ++order) {
            int nparam = alm_core->fcs->ndup[order].size();


            for (std::vector<ConstraintClass>::iterator p = alm_core->constraint->const_symmetry[order].begin();
                 p != alm_core->constraint->const_symmetry[order].end();
                 ++p) {
                ofs_fcs << "   0 = " << std::scientific << std::setprecision(6);
                ConstraintClass const_pointer = *p;
                for (j = 0; j < nparam; ++j) {
                    if (std::abs(const_pointer.w_const[j]) > eps8) {
                        str_tmp = " * (FC" + boost::lexical_cast<std::string>(order + 2)
                            + "_" + boost::lexical_cast<std::string>(j + 1) + ")";
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
    }

    ofs_fcs.unsetf(std::ios::showpos);

    for (order = 0; order < maxorder; ++order) {
        str_fcs[order] = "**FC" + boost::lexical_cast<std::string>(order + 2);
    }

    ofs_fcs << std::endl << std::endl;
    ofs_fcs << " ------------------------ All FCs below ------------------------" << std::endl;

    int ip = 0;
    int id;

    for (order = 0; order < maxorder; ++order) {

        id = 0;

        if (alm_core->fcs->ndup[order].size() > 0) {
            ofs_fcs << std::endl << std::setw(6) << str_fcs[order] << std::endl;

            for (unsigned int iuniq = 0; iuniq < alm_core->fcs->ndup[order].size(); ++iuniq) {

                str_tmp = "  # FC" + boost::lexical_cast<std::string>(order + 2) + "_";
                str_tmp += boost::lexical_cast<std::string>(iuniq + 1);

                ofs_fcs << str_tmp << std::setw(5) << alm_core->fcs->ndup[order][iuniq]
                    << std::setw(16) << std::scientific
                    << std::setprecision(7) << alm_core->fitting->params[ip] << std::endl;

                for (j = 0; j < alm_core->fcs->ndup[order][iuniq]; ++j) {
                    ofs_fcs << std::setw(5) << j + 1 << std::setw(12)
                        << std::setprecision(5) << std::fixed << alm_core->fcs->fc_table[order][id].coef;
                    for (k = 0; k < order + 2; ++k) {
                        ofs_fcs << std::setw(6)
                            << alm_core->fcs->easyvizint(alm_core->fcs->fc_table[order][id].elems[k]);
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

    std::cout << " Force constants in a human-readable format : "
        << alm_core->files->file_fcs << std::endl;
}

void Writer::write_displacement_pattern(ALM *alm)
{
    int i, j;
    int order;
    int counter;

    std::ofstream ofs_pattern;

    ALMCore *alm_core = alm->get_alm_core();
    int maxorder = alm_core->interaction->maxorder;

    std::cout << " Suggested displacement patterns are printed in the following files: " << std::endl;

    for (order = 0; order < maxorder; ++order) {
        ofs_pattern.open(alm_core->files->file_disp_pattern[order].c_str(), std::ios::out);
        if (!ofs_pattern)
            alm_core->error->exit("write_displacement_pattern",
                                  "Cannot open file_disp_pattern");

        counter = 0;

        ofs_pattern << "Basis : " << alm_core->displace->disp_basis[0] << std::endl;

        for (auto it = alm_core->displace->pattern_all[order].begin();
             it != alm_core->displace->pattern_all[order].end(); ++it) {
            AtomWithDirection entry = *it;

            ++counter;

            ofs_pattern << std::setw(5) << counter << ":"
                << std::setw(5) << entry.atoms.size() << std::endl;
            for (i = 0; i < entry.atoms.size(); ++i) {
                ofs_pattern << std::setw(7) << entry.atoms[i] + 1;
                for (j = 0; j < 3; ++j) {
                    ofs_pattern << std::setw(15) << entry.directions[3 * i + j];
                }
                ofs_pattern << std::endl;
            }
        }

        ofs_pattern.close();

        std::cout << "  " << alm_core->interaction->str_order[order]
            << " : " << alm_core->files->file_disp_pattern[order] << std::endl;
    }
    std::cout << std::endl;
}

void Writer::write_misc_xml(ALM *alm)
{
    SystemInfo system_structure;

    int i, j;

    ALMCore *alm_core = alm->get_alm_core();

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            system_structure.lattice_vector[i][j] = alm_core->system->lavec[i][j];
        }
    }

    system_structure.nat = alm_core->system->nat;
    system_structure.natmin = alm_core->symmetry->nat_prim;
    system_structure.ntran = alm_core->symmetry->ntran;
    system_structure.nspecies = alm_core->system->nkd;

    AtomProperty prop_tmp;

    for (i = 0; i < alm_core->system->nat; ++i) {
        prop_tmp.x = alm_core->system->xcoord[i][0];
        prop_tmp.y = alm_core->system->xcoord[i][1];
        prop_tmp.z = alm_core->system->xcoord[i][2];
        prop_tmp.kind = alm_core->system->kd[i];
        prop_tmp.atom = alm_core->symmetry->map_s2p[i].atom_num + 1;
        prop_tmp.tran = alm_core->symmetry->map_s2p[i].tran_num + 1;

        system_structure.atoms.push_back(AtomProperty(prop_tmp));
    }

    using boost::property_tree::ptree;

    ptree pt;
    std::string str_pos[3];

    pt.put("Data.ALM_version", ALAMODE_VERSION);
    pt.put("Data.Fitting.DisplaceFile", alm_core->files->file_disp);
    pt.put("Data.Fitting.ForceFile", alm_core->files->file_force);
    pt.put("Data.Fitting.Constraint", alm_core->constraint->constraint_mode);

    pt.put("Data.Structure.NumberOfAtoms", system_structure.nat);
    pt.put("Data.Structure.NumberOfElements", system_structure.nspecies);

    for (i = 0; i < system_structure.nspecies; ++i) {
        ptree &child = pt.add("Data.Structure.AtomicElements.element",
                              alm_core->system->kdname[i]);
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
    ss << alm_core->interaction->is_periodic[0] << " "
        << alm_core->interaction->is_periodic[1] << " "
        << alm_core->interaction->is_periodic[2];
    pt.put("Data.Structure.Periodicity", ss.str());

    pt.put("Data.Structure.Position", "");
    std::string str_tmp;

    for (i = 0; i < system_structure.nat; ++i) {
        str_tmp.clear();
        for (j = 0; j < 3; ++j) str_tmp += " " + double2string(alm_core->system->xcoord[i][j]);
        ptree &child = pt.add("Data.Structure.Position.pos", str_tmp);
        child.put("<xmlattr>.index", i + 1);
        child.put("<xmlattr>.element", alm_core->system->kdname[alm_core->system->kd[i] - 1]);
    }

    pt.put("Data.Symmetry.NumberOfTranslations", alm_core->symmetry->ntran);
    for (i = 0; i < system_structure.ntran; ++i) {
        for (j = 0; j < system_structure.natmin; ++j) {
            ptree &child = pt.add("Data.Symmetry.Translations.map",
                                  alm_core->symmetry->map_p2s[j][i] + 1);
            child.put("<xmlattr>.tran", i + 1);
            child.put("<xmlattr>.atom", j + 1);
        }
    }

    if (alm_core->system->lspin) {
        pt.put("Data.MagneticMoments", "");
        pt.put("Data.MagneticMoments.Noncollinear", alm_core->system->noncollinear);
        pt.put("Data.MagneticMoments.TimeReversalSymmetry", alm_core->symmetry->trev_sym_mag);
        for (i = 0; i < system_structure.nat; ++i) {
            str_tmp.clear();
            for (j = 0; j < 3; ++j) str_tmp += " " + double2string(alm_core->system->magmom[i][j], 5);
            ptree &child = pt.add("Data.MagneticMoments.mag", str_tmp);
            child.put("<xmlattr>.index", i + 1);
        }
    }

    pt.put("Data.ForceConstants", "");
    str_tmp.clear();

    pt.put("Data.ForceConstants.HarmonicUnique.NFC2", alm_core->fcs->ndup[0].size());

    int ihead = 0;
    int k = 0;
    int nelem = alm_core->interaction->maxorder + 1;
    int *pair_tmp;

    allocate(pair_tmp, nelem);

    for (unsigned int ui = 0; ui < alm_core->fcs->ndup[0].size(); ++ui) {

        for (i = 0; i < 2; ++i) {
            pair_tmp[i] = alm_core->fcs->fc_table[0][ihead].elems[i] / 3;
        }
        j = alm_core->symmetry->map_s2p[pair_tmp[0]].atom_num;

        ptree &child = pt.add("Data.ForceConstants.HarmonicUnique.FC2",
                              double2string(alm_core->fitting->params[k]));
        child.put("<xmlattr>.pairs",
                  boost::lexical_cast<std::string>(alm_core->fcs->fc_table[0][ihead].elems[0])
                  + " " + boost::lexical_cast<std::string>(alm_core->fcs->fc_table[0][ihead].elems[1]));
        child.put("<xmlattr>.multiplicity",
                  alm_core->interaction->mindist_pairs[pair_tmp[0]][pair_tmp[1]].size());
        ihead += alm_core->fcs->ndup[0][ui];
        ++k;
    }
    ihead = 0;

    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<MinimumDistanceCluster>::iterator iter_cluster;
    int multiplicity;

    if (alm_core->interaction->maxorder > 1) {

        pt.put("Data.ForceConstants.CubicUnique.NFC3", alm_core->fcs->ndup[1].size());

        for (unsigned int ui = 0; ui < alm_core->fcs->ndup[1].size(); ++ui) {
            for (i = 0; i < 3; ++i) {
                pair_tmp[i] = alm_core->fcs->fc_table[1][ihead].elems[i] / 3;
            }
            j = alm_core->symmetry->map_s2p[pair_tmp[0]].atom_num;

            atom_tmp.clear();
            for (i = 1; i < 3; ++i) {
                atom_tmp.push_back(pair_tmp[i]);
            }
            std::sort(atom_tmp.begin(), atom_tmp.end());

            iter_cluster = alm_core->interaction->mindist_cluster[1][j].find(
                MinimumDistanceCluster(atom_tmp, cell_dummy));
            if (iter_cluster == alm_core->interaction->mindist_cluster[1][j].end()) {
                alm_core->error->exit("load_reference_system_xml",
                                      "Cubic force constant is not found.");
            } else {
                multiplicity = (*iter_cluster).cell.size();
            }

            ptree &child = pt.add("Data.ForceConstants.CubicUnique.FC3",
                                  double2string(alm_core->fitting->params[k]));
            child.put("<xmlattr>.pairs",
                      boost::lexical_cast<std::string>(alm_core->fcs->fc_table[1][ihead].elems[0])
                      + " " + boost::lexical_cast<std::string>(alm_core->fcs->fc_table[1][ihead].elems[1])
                      + " " + boost::lexical_cast<std::string>(alm_core->fcs->fc_table[1][ihead].elems[2]));
            child.put("<xmlattr>.multiplicity", multiplicity);
            ihead += alm_core->fcs->ndup[1][ui];
            ++k;
        }
    }

    int ip, ishift;

    std::sort(alm_core->fcs->fc_table[0].begin(), alm_core->fcs->fc_table[0].end());

    for (auto it = alm_core->fcs->fc_table[0].begin();
         it != alm_core->fcs->fc_table[0].end(); ++it) {
        FcProperty fctmp = *it;
        ip = fctmp.mother;

        for (k = 0; k < 2; ++k) {
            pair_tmp[k] = fctmp.elems[k] / 3;
        }
        j = alm_core->symmetry->map_s2p[pair_tmp[0]].atom_num;
        for (auto it2 = alm_core->interaction->mindist_pairs[pair_tmp[0]][pair_tmp[1]].begin();
             it2 != alm_core->interaction->mindist_pairs[pair_tmp[0]][pair_tmp[1]].end(); ++it2) {
            ptree &child = pt.add("Data.ForceConstants.HARMONIC.FC2",
                                  double2string(alm_core->fitting->params[ip] * fctmp.coef
                                      / static_cast<double>(alm_core->interaction->mindist_pairs[pair_tmp[0]][pair_tmp[1]].size())));

            child.put("<xmlattr>.pair1", boost::lexical_cast<std::string>(j + 1)
                      + " " + boost::lexical_cast<std::string>(fctmp.elems[0] % 3 + 1));
            child.put("<xmlattr>.pair2", boost::lexical_cast<std::string>(pair_tmp[1] + 1)
                      + " " + boost::lexical_cast<std::string>(fctmp.elems[1] % 3 + 1)
                      + " " + boost::lexical_cast<std::string>((*it2).cell + 1));
        }
    }

    ishift = alm_core->fcs->ndup[0].size();

    // Print anharmonic force constants to the xml file.

    int imult;

    int order;
    std::string elementname;
    for (order = 1; order < alm_core->interaction->maxorder; ++order) {

        std::sort(alm_core->fcs->fc_table[order].begin(), alm_core->fcs->fc_table[order].end());

        for (auto it = alm_core->fcs->fc_table[order].begin();
             it != alm_core->fcs->fc_table[order].end(); ++it) {
            FcProperty fctmp = *it;
            ip = fctmp.mother + ishift;

            for (k = 0; k < order + 2; ++k) {
                pair_tmp[k] = fctmp.elems[k] / 3;
            }
            j = alm_core->symmetry->map_s2p[pair_tmp[0]].atom_num;

            atom_tmp.clear();

            for (k = 1; k < order + 2; ++k) {
                atom_tmp.push_back(pair_tmp[k]);
            }
            std::sort(atom_tmp.begin(), atom_tmp.end());

            elementname = "Data.ForceConstants.ANHARM"
                + boost::lexical_cast<std::string>(order + 2)
                + ".FC" + boost::lexical_cast<std::string>(order + 2);


            iter_cluster = alm_core->interaction->mindist_cluster[order][j].find(
                MinimumDistanceCluster(atom_tmp, cell_dummy));

            if (iter_cluster != alm_core->interaction->mindist_cluster[order][j].end()) {
                multiplicity = (*iter_cluster).cell.size();

                for (imult = 0; imult < multiplicity; ++imult) {
                    std::vector<int> cell_now = (*iter_cluster).cell[imult];

                    ptree &child = pt.add(elementname,
                                          double2string(alm_core->fitting->params[ip] * fctmp.coef
                                              / static_cast<double>(multiplicity)));

                    child.put("<xmlattr>.pair1", boost::lexical_cast<std::string>(j + 1)
                              + " " + boost::lexical_cast<std::string>(fctmp.elems[0] % 3 + 1));

                    for (k = 1; k < order + 2; ++k) {
                        child.put("<xmlattr>.pair" + boost::lexical_cast<std::string>(k + 1),
                                  boost::lexical_cast<std::string>(pair_tmp[k] + 1)
                                  + " " + boost::lexical_cast<std::string>(fctmp.elems[k] % 3 + 1)
                                  + " " + boost::lexical_cast<std::string>(cell_now[k - 1] + 1));
                    }
                }
            } else {
                alm_core->error->exit("write_misc_xml", "This cannot happen.");
            }
        }
        ishift += alm_core->fcs->ndup[order].size();
    }

    using namespace boost::property_tree::xml_parser;
    const int indent = 2;

    std::string file_xml = alm_core->files->job_title + ".xml";

#if BOOST_VERSION >= 105600
    write_xml(file_xml, pt, std::locale(),
              xml_writer_make_settings<ptree::key_type>(' ', indent,
                                                        widen<std::string>("utf-8")));
#else
    write_xml(file_xml, pt, std::locale(),
              xml_writer_make_settings(' ', indent, widen<char>("utf-8")));
#endif

    deallocate(pair_tmp);

    std::cout << " Input data for the phonon code ANPHON      : " << file_xml << std::endl;
}

void Writer::write_hessian(ALM *alm)
{
    int i, j, itran, ip;
    int pair_tmp[2];
    int pair_tran[2];
    std::ofstream ofs_hes;
    double **hessian;

    ALMCore *alm_core = alm->get_alm_core();
    int nat3 = 3 * alm_core->system->nat;

    allocate(hessian, nat3, nat3);

    for (i = 0; i < nat3; ++i) {
        for (j = 0; j < nat3; ++j) {
            hessian[i][j] = 0.0;
        }
    }

    for (auto it = alm_core->fcs->fc_table[0].begin();
         it != alm_core->fcs->fc_table[0].end(); ++it) {
        FcProperty fctmp = *it;
        ip = fctmp.mother;

        for (i = 0; i < 2; ++i) pair_tmp[i] = fctmp.elems[i] / 3;
        for (itran = 0; itran < alm_core->symmetry->ntran; ++itran) {
            for (i = 0; i < 2; ++i) {
                pair_tran[i] = alm_core->symmetry->map_sym[pair_tmp[i]][alm_core->symmetry->symnum_tran[itran]];
            }
            hessian[3 * pair_tran[0] + fctmp.elems[0] % 3][3 * pair_tran[1] + fctmp.elems[1] % 3]
                = alm_core->fitting->params[ip] * fctmp.coef;
        }
    }

    ofs_hes.open(alm_core->files->file_hes.c_str(), std::ios::out);
    if (!ofs_hes) alm_core->error->exit("write_hessian", "cannot create hessian file");

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

    std::cout << " Complete Hessian matrix                    : " << alm_core->files->file_hes << std::endl;
}

std::string Writer::double2string(const double d, const int nprec)
{
    std::string rt;
    std::stringstream ss;

    ss << std::scientific << std::setprecision(nprec) << d;
    ss >> rt;
    return rt;
}

void Writer::write_in_QEformat(ALMCore *alm)
{
    int i, j, itran, ip;
    int pair_tmp[2];
    int pair_tran[2];
    std::ofstream ofs_hes;
    double **hessian;
    int nat3 = 3 * alm->system->nat;

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
                = alm->fitting->params[ip] * fctmp.coef;
        }
    }

    std::string file_fc = alm->files->job_title + ".fc";

    ofs_hes.open(file_fc.c_str(), std::ios::out);
    if (!ofs_hes) alm->error->exit("write_hessian", "cannot create hessian file");

    ofs_hes << "  1  1  1" << std::endl;
    for (int icrd = 0; icrd < 3; ++icrd) {
        for (int jcrd = 0; jcrd < 3; ++jcrd) {
            for (i = 0; i < alm->system->nat; ++i) {
                for (j = 0; j < alm->system->nat; ++j) {
                    ofs_hes << std::setw(3) << icrd + 1;
                    ofs_hes << std::setw(3) << jcrd + 1;
                    ofs_hes << std::setw(3) << i + 1;
                    ofs_hes << std::setw(3) << j + 1;
                    ofs_hes << std::endl;
                    ofs_hes << "  1  1  1 " << std::setw(20) << std::setprecision(13) <<  std::scientific << hessian[3 * j + jcrd][3 * i + icrd];
                    ofs_hes << std::endl;
                }
            }
        }
    }
    ofs_hes.close();
    deallocate(hessian);
}
