/*
 system.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "system.h"
#include "constants.h"
#include "constraint.h"
#include "error.h"
#include "fcs.h"
#include "fitting.h"
#include "mathfunctions.h"
#include "memory.h"
#include "symmetry.h"
#include "timer.h"
#include "xml_parser.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

System::System()
{
    set_default_variables();
}

System::~System()
{
    deallocate_variables();
}

void System::init(ALM *alm)
{
    int i, j;
    alm->timer->start_clock("system");

    set_cell(lavec, nat, nkd, kd, xcoord, supercell);
    setup_atomic_class(kd);

    print_structure_stdout(supercell);
    if (lspin) print_magmom_stdout();
    alm->timer->print_elapsed();
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    alm->timer->stop_clock("system");
}

void System::set_cell(const double lavec_in[3][3],
                      const unsigned int nat_in,
                      const unsigned int nkd_in,
                      int *kd_in,
                      double **xf_in,
                      Cell &cell_out)
{
    unsigned int i, j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            cell_out.lattice_vector[i][j] = lavec_in[i][j];
        }
    }
    set_reciprocal_latt(cell_out.lattice_vector,
                        cell_out.reciprocal_lattice_vector);

    cell_out.volume = volume(cell_out.lattice_vector, Direct);
    cell_out.number_of_atmos = nat_in;
    cell_out.number_of_elems = nkd_in;
    cell_out.kind.clear();
    cell_out.x_fractional.clear();
    cell_out.x_cartesian.clear();

    std::vector<double> xtmp;
    xtmp.resize(3);
    for (i = 0; i < nat_in; ++i) {
        cell_out.kind.push_back(kd_in[i]);
        for (j = 0; j < 3; ++j) {
            xtmp[j] = xf_in[i][j];
        }
        cell_out.x_fractional.push_back(xtmp);
    }

    double xf_tmp[3], xc_tmp[3];

    for (const auto &xf : cell_out.x_fractional) {
        for (i = 0; i < 3; ++i) {
            xf_tmp[i] = xf[i];
        }
        rotvec(xc_tmp, xf_tmp, cell_out.lattice_vector);
        for (i = 0; i < 3; ++i) {
            xtmp[i] = xc_tmp[i];
        }
        cell_out.x_cartesian.push_back(xtmp);
    }
}

void System::set_reciprocal_latt(const double aa[3][3], double bb[3][3])
{
    /*
    Calculate Reciprocal Lattice Vectors

    Here, BB is just the inverse matrix of AA (multiplied by factor 2 Pi)

    BB = 2 Pi AA^{-1},
    = t(b1, b2, b3)

    (b11 b12 b13)
    = (b21 b22 b23)
    (b31 b32 b33),

    b1 = t(b11, b12, b13) etc.
    */

    double det;
    det = aa[0][0] * aa[1][1] * aa[2][2]
        + aa[1][0] * aa[2][1] * aa[0][2]
        + aa[2][0] * aa[0][1] * aa[1][2]
        - aa[0][0] * aa[2][1] * aa[1][2]
        - aa[2][0] * aa[1][1] * aa[0][2]
        - aa[1][0] * aa[0][1] * aa[2][2];

    if (std::abs(det) < eps12) {
        exit("set_reciprocal_latt", "Lattice Vector is singular");
    }

    double factor = 2.0 * pi / det;

    bb[0][0] = (aa[1][1] * aa[2][2] - aa[1][2] * aa[2][1]) * factor;
    bb[0][1] = (aa[0][2] * aa[2][1] - aa[0][1] * aa[2][2]) * factor;
    bb[0][2] = (aa[0][1] * aa[1][2] - aa[0][2] * aa[1][1]) * factor;

    bb[1][0] = (aa[1][2] * aa[2][0] - aa[1][0] * aa[2][2]) * factor;
    bb[1][1] = (aa[0][0] * aa[2][2] - aa[0][2] * aa[2][0]) * factor;
    bb[1][2] = (aa[0][2] * aa[1][0] - aa[0][0] * aa[1][2]) * factor;

    bb[2][0] = (aa[1][0] * aa[2][1] - aa[1][1] * aa[2][0]) * factor;
    bb[2][1] = (aa[0][1] * aa[2][0] - aa[0][0] * aa[2][1]) * factor;
    bb[2][2] = (aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]) * factor;
}

void System::frac2cart(double **xf)
{
    // x_cartesian = A x_fractional

    int i, j;

    double *x_tmp;
    allocate(x_tmp, 3);

    for (i = 0; i < nat; ++i) {

        rotvec(x_tmp, xf[i], lavec);

        for (j = 0; j < 3; ++j) {
            xf[i][j] = x_tmp[j];
        }
    }
    deallocate(x_tmp);
}

void System::load_reference_system_xml(Symmetry *symmetry,
                                       Fcs *fcs,
                                       std::string file_reference_fcs,
                                       const int order_fcs,
                                       double *const_out)
{
    using namespace boost::property_tree;
    ptree pt;

    int nat_ref, natmin_ref, ntran_ref;
    int **intpair_ref;
    std::string str_error;
    double *fcs_ref;
    int nfcs_ref;

    try {
        read_xml(file_reference_fcs, pt);
    }
    catch (std::exception &e) {
        if (order_fcs == 0) {
            str_error = "Cannot open file FC2XML ( " + file_reference_fcs + " )";
        } else if (order_fcs == 1) {
            str_error = "Cannot open file FC3XML ( " + file_reference_fcs + " )";
        }
        exit("load_reference_system_xml", str_error.c_str());
    }

    nat_ref = boost::lexical_cast<int>(
        get_value_from_xml(pt, "Data.Structure.NumberOfAtoms"));
    ntran_ref = boost::lexical_cast<int>(
        get_value_from_xml(pt, "Data.Symmetry.NumberOfTranslations"));
    natmin_ref = nat_ref / ntran_ref;

    if (natmin_ref != symmetry->nat_prim) {
        exit("load_reference_system_xml",
             "The number of atoms in the primitive cell is not consistent.");
    }

    if (order_fcs == 0) {
        nfcs_ref = boost::lexical_cast<int>(
            get_value_from_xml(pt, "Data.ForceConstants.HarmonicUnique.NFC2"));

        if (nfcs_ref != fcs->nequiv[0].size()) {
            exit("load_reference_system_xml",
                 "The number of harmonic force constants is not the same.");
        }

    } else if (order_fcs == 1) {
        nfcs_ref = boost::lexical_cast<int>(
            get_value_from_xml(pt, "Data.ForceConstants.CubicUnique.NFC3"));

        if (nfcs_ref != fcs->nequiv[1].size()) {
            exit("load_reference_system_xml",
                 "The number of cubic force constants is not the same.");
        }
    }
    allocate(fcs_ref, nfcs_ref);
    allocate(intpair_ref, nfcs_ref, 3);

    int counter = 0;

    if (order_fcs == 0) {
        BOOST_FOREACH (const ptree::value_type& child_, pt.get_child("Data.ForceConstants.HarmonicUnique")) {
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
    } else if (order_fcs == 1) {
        BOOST_FOREACH (const ptree::value_type& child_, pt.get_child("Data.ForceConstants.CubicUnique")) {
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

    int i;
    std::set<FcProperty> list_found;
    std::set<FcProperty>::iterator iter_found;
    int *ind;
    int nterms = order_fcs + 2;
    allocate(ind, nterms);

    list_found.clear();

    for (auto p = fcs->fc_table[order_fcs].begin(); p != fcs->fc_table[order_fcs].end(); ++p) {
        FcProperty list_tmp = *p; // Using copy constructor
        for (i = 0; i < nterms; ++i) {
            ind[i] = list_tmp.elems[i];
        }
        list_found.insert(FcProperty(nterms, list_tmp.sign,
                                     ind, list_tmp.mother));
    }

    for (i = 0; i < nfcs_ref; ++i) {
        iter_found = list_found.find(FcProperty(nterms, 1.0,
                                                intpair_ref[i], 1));
        if (iter_found == list_found.end()) {
            exit("load_reference_system",
                 "Cannot find equivalent force constant, number: ",
                 i + 1);
        }
        FcProperty arrtmp = *iter_found;
        const_out[arrtmp.mother] = fcs_ref[i];
    }

    deallocate(intpair_ref);
    deallocate(fcs_ref);
    deallocate(ind);
    list_found.clear();
}


double System::volume(const double latt_in[3][3], LatticeType type)
{
    int i, j;
    double mat[3][3];

    if (type == Direct) {
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                mat[i][j] = latt_in[j][i];
            }
        }
    } else if (type == Reciprocal) {
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                mat[i][j] = latt_in[i][j];
            }
        }
    } else {
        exit("volume", "Invalid LatticeType is given");
    }

    auto vol = std::abs(mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
        + mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2])
        + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]));

    return vol;
}


void System::set_default_variables()
{
    // Default values
    nat = 0;
    nkd = 0;
    kd = nullptr;
    kdname = nullptr;
    ndata = 0;
    nstart = 1;
    nend = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lavec[i][j] = 0;
        }
    }
    xcoord = nullptr;
    str_magmom = "";
    atomlist_class = nullptr;
    nclassatom = 0;
    lspin = false;
    noncollinear = 0;
    magmom = nullptr;
    trev_sym_mag = 1;
}

void System::deallocate_variables()
{
    if (kd) {
        deallocate(kd);
    }
    if (kdname) {
        deallocate(kdname);
    }
    if (xcoord) {
        deallocate(xcoord);
    }
    if (atomlist_class) {
        deallocate(atomlist_class);
    }
    if (magmom) {
        deallocate(magmom);
    }
}

void System::setup_atomic_class(int *kd)
{
    // In the case of collinear calculation, spin moments are considered as scalar
    // variables. Therefore, the same elements with different magnetic moments are
    // considered as different types. In noncollinear calculations, 
    // magnetic moments are not considered in this stage. They will be treated
    // separately in symmetry.cpp where spin moments will be rotated and flipped 
    // using time-reversal symmetry.

    unsigned int i;
    AtomType type_tmp;
    std::set<AtomType> set_type;
    set_type.clear();

    for (i = 0; i < nat; ++i) {
        type_tmp.element = kd[i];

        if (noncollinear == 0) {
            type_tmp.magmom = magmom[i][2];
        } else {
            type_tmp.magmom = 0.0;
        }
        set_type.insert(type_tmp);
    }

    nclassatom = set_type.size();

    allocate(atomlist_class, nclassatom);

    for (i = 0; i < nat; ++i) {
        int count = 0;
        for (auto it = set_type.cbegin(); it != set_type.cend(); ++it) {
            if (noncollinear) {
                if (kd[i] == (*it).element) {
                    atomlist_class[count].push_back(i);
                }
            } else {
                if ((kd[i] == (*it).element)
                    && (std::abs(magmom[i][2] - (*it).magmom) < eps6)) {
                    atomlist_class[count].push_back(i);
                }
            }
            ++count;
        }
    }
    set_type.clear();
}


void System::print_structure_stdout(const Cell &cell)
{
    using namespace std;
    int i, j;

    cout << " SYSTEM" << endl;
    cout << " ======" << endl << endl;

    cout.setf(ios::scientific);

    cout << "  Lattice Vector" << endl;
    cout << setw(16) << cell.lattice_vector[0][0];
    cout << setw(15) << cell.lattice_vector[1][0];
    cout << setw(15) << cell.lattice_vector[2][0];
    cout << " : a1" << endl;

    cout << setw(16) << cell.lattice_vector[0][1];
    cout << setw(15) << cell.lattice_vector[1][1];
    cout << setw(15) << cell.lattice_vector[2][1];
    cout << " : a2" << endl;

    cout << setw(16) << cell.lattice_vector[0][2];
    cout << setw(15) << cell.lattice_vector[1][2];
    cout << setw(15) << cell.lattice_vector[2][2];
    cout << " : a3" << endl;
    cout << endl;

    cout << "  Cell volume = " << cell.volume << " (a.u)^3"
        << endl << endl;

    cout << "  Reciprocal Lattice Vector" << std::endl;
    cout << setw(16) << supercell.reciprocal_lattice_vector[0][0];
    cout << setw(15) << supercell.reciprocal_lattice_vector[0][1];
    cout << setw(15) << supercell.reciprocal_lattice_vector[0][2];
    cout << " : b1" << endl;

    cout << setw(16) << supercell.reciprocal_lattice_vector[1][0];
    cout << setw(15) << supercell.reciprocal_lattice_vector[1][1];
    cout << setw(15) << supercell.reciprocal_lattice_vector[1][2];
    cout << " : b2" << endl;

    cout << setw(16) << supercell.reciprocal_lattice_vector[2][0];
    cout << setw(15) << supercell.reciprocal_lattice_vector[2][1];
    cout << setw(15) << supercell.reciprocal_lattice_vector[2][2];
    cout << " : b3" << endl;
    cout << endl;

    cout << "  Atomic species:" << endl;
    for (i = 0; i < cell.number_of_elems; ++i) {
        cout << setw(6) << i + 1 << setw(5) << kdname[i] << endl;
    }
    cout << endl;

    cout << "  Atomic positions in fractional basis and atomic species" << endl;
    for (i = 0; i < cell.number_of_atmos; ++i) {
        cout << setw(6) << i + 1;
        cout << setw(15) << cell.x_fractional[i][0];
        cout << setw(15) << cell.x_fractional[i][1];
        cout << setw(15) << cell.x_fractional[i][2];
        cout << setw(5) << cell.kind[i] << endl;
    }
    cout << endl << endl;
    cout.unsetf(ios::scientific);
}


void System::print_magmom_stdout()
{
    using namespace std;

    cout << "  MAGMOM is given. The magnetic moments of each atom are as follows:" << endl;
    for (auto i = 0; i < supercell.number_of_atmos; ++i) {
        cout << setw(6) << i + 1;
        cout << setw(5) << magmom[i][0];
        cout << setw(5) << magmom[i][1];
        cout << setw(5) << magmom[i][2];
        cout << endl;
    }
    cout << endl;
    if (noncollinear == 0) {
        cout << "  NONCOLLINEAR = 0: magnetic moments are considered as scalar variables." << endl;
    } else if (noncollinear == 1) {
        cout << "  NONCOLLINEAR = 1: magnetic moments are considered as vector variables." << endl;
        if (trev_sym_mag) {
            cout << "  TREVSYM = 1: Time-reversal symmetry will be considered for generating magnetic space group"
                << endl;
        } else {
            cout <<
                "  TREVSYM = 0: Time-reversal symmetry will NOT be considered for generating magnetic space group"
                << endl;
        }
    }
    cout << endl << endl;
}
