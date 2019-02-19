/*
 system.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "system.h"
#include "constants.h"
#include "error.h"
#include "mathfunctions.h"
#include "memory.h"
#include "timer.h"
#include <iostream>
#include <iomanip>
#include <set>

using namespace ALM_NS;

System::System()
{
    set_default_variables();
}

System::~System()
{
    deallocate_variables();
}

void System::init(const int verbosity,
                  Timer *timer)
{
    const auto nat = supercell.number_of_atoms;

    timer->start_clock("system");

    // Set atomic types (kind + magmom)
    set_atomtype_group();

    const auto nneib = 27;
    if (x_image) {
        deallocate(x_image);
    }
    allocate(x_image, nneib, nat, 3);

    if (exist_image) {
        deallocate(exist_image);
    }
    allocate(exist_image, nneib);

    generate_coordinate_of_periodic_images();

    if (verbosity > 0) {
        print_structure_stdout(supercell);
        if (spin.lspin) print_magmom_stdout();
        timer->print_elapsed();
        std::cout << " -------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
    }

    timer->stop_clock("system");
}

void System::set_supercell(const double lavec_in[3][3],
                           const size_t nat_in,
                           const int *kind_in,
                           const double xf_in[][3],
                           const double transformation_matrix[3][3],
                           const double primitive_axes[3][3])
{
    set_inputcell(lavec_in, nat_in, kind_in, xf_in);
    build_supercell(transformation_matrix);
}

void System::set_inputcell(const double lavec_in[3][3],
                           const size_t nat_in,
                           const int *kind_in,
                           const double xf_in[][3])
{
    size_t i, j;
    std::vector<int> unique_nums(nat_in);
    auto wrong_number = false;
    bool in_unique_nums;

    for (i = 0; i < nat_in; i++) {
        unique_nums[i] = 0;
    }

    size_t nkd = 0;
    for (i = 0; i < nat_in; i++) {
        in_unique_nums = false;
        for (j = 0; j < nkd; j++) {
            if (unique_nums[j] == kind_in[i]) {
                in_unique_nums = true;
                break;
            }
        }
        if (!in_unique_nums) {
            unique_nums[nkd] = kind_in[i];
            nkd++;
        }
    }

    for (i = 0; i < nkd; i++) {
        if (static_cast<size_t>(unique_nums[i]) > nkd) {
            std::cout << " WARNING : integers assigned to atoms are wrong. "
                << " The numbers will be resorted." << std::endl;
            wrong_number = true;
            break;
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            inputcell.lattice_vector[i][j] = lavec_in[i][j];
        }
    }
    set_reciprocal_latt(inputcell.lattice_vector,
                        inputcell.reciprocal_lattice_vector);

    inputcell.volume = volume(inputcell.lattice_vector, Direct);
    inputcell.number_of_atoms = nat_in;
    inputcell.number_of_elems = nkd;
    inputcell.kind.clear();
    inputcell.kind.shrink_to_fit();
    inputcell.x_fractional.clear();
    inputcell.x_fractional.shrink_to_fit();
    inputcell.x_cartesian.clear();
    inputcell.x_cartesian.shrink_to_fit();

    std::vector<double> xtmp;

    if (!wrong_number) {
        for (i = 0; i < nat_in; ++i) {
            inputcell.kind.push_back(kind_in[i]);
        }
    } else {
        for (i = 0; i < nat_in; ++i) {
            for (j = 0; j < nkd; j++) {
                if (kind_in[i] == unique_nums[j]) {
                    inputcell.kind.push_back(static_cast<int>(j + 1));
                }
            }
        }
    }

    xtmp.resize(3);
    for (i = 0; i < nat_in; ++i) {
        for (j = 0; j < 3; ++j) {
            xtmp[j] = xf_in[i][j];
        }
        inputcell.x_fractional.push_back(xtmp);
    }

    for (const auto &xf : inputcell.x_fractional) {
        rotvec(&xtmp[0], &xf[0], inputcell.lattice_vector);
        inputcell.x_cartesian.push_back(xtmp);
    }

    // This is needed to avoid segmentation fault.
    inputspin.magmom.clear();
    std::vector<double> vec(3);
    for (i = 0; i < nat_in; ++i) {
        for (j = 0; j < 3; ++j) {
            vec[j] = 0;
        }
        inputspin.magmom.push_back(vec);
    }
}

void System::build_supercell(const double transformation_matrix[3][3])
{
    // Build supercell from inputcell and transformation matrix.
    // Also, spin related variables of supercell is generated from inputspin.

    // Convention of the transformation matrix (Transmat)
    // (a_s, b_s, c_s) = (a_in, b_in, c_in) * Transmat

    matmul3(supercell.lattice_vector, 
            inputcell.lattice_vector, 
            transformation_matrix);

    set_reciprocal_latt(supercell.lattice_vector,
                        supercell.reciprocal_lattice_vector);

    supercell.volume = volume(supercell.lattice_vector, Direct);

    double tmat_determinant = determinant3(transformation_matrix);
    auto ncells = nint(std::abs(tmat_determinant));

    supercell.number_of_atoms = inputcell.number_of_atoms * ncells;
    supercell.number_of_elems = inputcell.number_of_elems;
    supercell.kind.clear();
    supercell.kind.shrink_to_fit();
    supercell.x_fractional.clear();
    supercell.x_fractional.shrink_to_fit();
    supercell.x_cartesian.clear();
    supercell.x_cartesian.shrink_to_fit();

    double convmat[3][3];

    matmul3(convmat, 
            supercell.reciprocal_lattice_vector, 
            inputcell.lattice_vector);

    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            convmat[i][j] /= 2.0 * pi;
        }
    }

    int nsize[3];

    for (auto i = 0; i < 3; ++i) {
        nsize[i] = 0;
        auto nmin = 0;
        auto nmax = 0;
        for (auto j = 0; j < 3; ++j) {
            nmin = std::min<int>(nmin, nint(transformation_matrix[i][j]));
            nmax = std::max<int>(nmax, nint(transformation_matrix[i][j]));
        }
        nsize[i] = nmax - nmin;
    }

    double xfrac[3], xshift[3];
    std::vector<double> xfrac_s(3);
    std::vector<std::vector<double>> x_shifted;

    bool new_entry, found_all;

    for (auto iat = 0; iat < inputcell.number_of_atoms; ++iat) {
        for (auto i = 0; i < 3; ++i) {
            xfrac[i] = inputcell.x_fractional[iat][i];
        }

        x_shifted.clear();
        x_shifted.shrink_to_fit();
        found_all = false;

        for (auto isize = 0; isize < nsize[0]; ++isize) {
            xshift[0] = xfrac[0] + static_cast<double>(isize);
            for (auto jsize = 0; jsize < nsize[1]; ++jsize) {
                xshift[1] = xfrac[1] + static_cast<double>(jsize);
                for (auto ksize = 0; ksize < nsize[2]; ++ksize) {
                    xshift[2] = xfrac[2] + static_cast<double>(ksize);
                    rotvec(&xfrac_s[0], xshift, convmat);

                    for (auto i = 0; i < 3; ++i) {
                        xfrac_s[i] -= static_cast<double>(nint(xfrac_s[i]));
                       while (xfrac_s[i] >= 1.0 - eps6) {
                           xfrac_s[i] -= 1.0;
                       }
                        while (xfrac_s[i] < -eps6) {
                            xfrac_s[i] += 1.0;
                        }
                    }
                    // std::cout << std::scientific;
                    // std::cout << std::setw(5) << iat + 1;
                    // std::cout << std::setw(15) << xfrac_s[0];
                    // std::cout << std::setw(15) << xfrac_s[1];
                    // std::cout << std::setw(15) << xfrac_s[2] << '\n';
                    
                    new_entry = true;
                    for (const auto &it: x_shifted) {
                        auto diff = std::pow(xfrac_s[0] - it[0], 2) 
                                  + std::pow(xfrac_s[1] - it[1], 2) 
                                  + std::pow(xfrac_s[2] - it[2], 2);

                        if (std::sqrt(diff) < eps8) {
                            new_entry = false;
                            break;
                        }
                    }

                    if (new_entry) x_shifted.push_back(xfrac_s);

                    if (x_shifted.size() == ncells) {
                        found_all = true;

                        for (const auto &it : x_shifted) {
                            supercell.x_fractional.push_back(it);
                            supercell.kind.push_back(inputcell.kind[iat]);
                            spin.magmom.push_back(inputspin.magmom[iat]);
                        }
                        break;
                    }
                }
                if (found_all) break;
            }
            if (found_all) break;
        }
    }

    std::vector<double> xtmp(3), xf(3);
    for (const auto &xf : supercell.x_fractional) {
        rotvec(&xtmp[0], &xf[0], supercell.lattice_vector);
        supercell.x_cartesian.push_back(xtmp);
    }
    spin.lspin = inputspin.lspin;
    spin.time_reversal_symm = inputspin.time_reversal_symm;
    spin.noncollinear = inputspin.noncollinear;
}

const Cell& System::get_supercell() const
{
    return supercell;
}

double*** System::get_x_image() const
{
    return x_image;
}

int* System::get_exist_image() const
{
    return exist_image;
}

void System::set_periodicity(const int is_periodic_in[3])
{
    if (! is_periodic) {
        // This should be already allocated though.
        allocate(is_periodic, 3);
    }
    for (unsigned int i = 0; i < 3; i++) {
        is_periodic[i] = is_periodic_in[i];
    }
}

int* System::get_periodicity() const
{
    return is_periodic;
}

void System::set_kdname(const std::string *kdname_in)
{
    const auto nkd = inputcell.number_of_elems;

    if (kdname) {
        deallocate(kdname);
    }
    allocate(kdname, nkd);
    for (size_t i = 0; i < nkd; ++i) {
        kdname[i] = kdname_in[i];
    }
}

std::string* System::get_kdname() const
{
    return kdname;
}

void System::set_reciprocal_latt(const double aa[3][3],
                                 double bb[3][3]) const
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

    const auto det
        = aa[0][0] * aa[1][1] * aa[2][2]
        + aa[1][0] * aa[2][1] * aa[0][2]
        + aa[2][0] * aa[0][1] * aa[1][2]
        - aa[0][0] * aa[2][1] * aa[1][2]
        - aa[2][0] * aa[1][1] * aa[0][2]
        - aa[1][0] * aa[0][1] * aa[2][2];

    if (std::abs(det) < eps12) {
        exit("set_reciprocal_latt", "Lattice Vector is singular");
    }

    const auto factor = 2.0 * pi / det;

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

void System::frac2cart(double **xf) const
{
    // x_cartesian = A x_fractional

    double *x_tmp;
    allocate(x_tmp, 3);

    for (size_t i = 0; i < supercell.number_of_atoms; ++i) {

        rotvec(x_tmp, xf[i], supercell.lattice_vector);

        for (auto j = 0; j < 3; ++j) {
            xf[i][j] = x_tmp[j];
        }
    }
    deallocate(x_tmp);
}

double System::volume(const double latt_in[3][3],
                      const LatticeType type) const
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

    const auto vol = std::abs(mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
        + mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2])
        + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]));

    return vol;
}

void System::set_default_variables()
{
    kdname = nullptr;

    inputcell.number_of_atoms = 0;
    inputcell.number_of_elems = 0;
    supercell.number_of_atoms = 0;
    supercell.number_of_elems = 0;
    primcell.number_of_atoms = 0;
    primcell.number_of_elems = 0;

    allocate(is_periodic, 3);
    is_periodic[0] = 1;
    is_periodic[1] = 1;
    is_periodic[2] = 1;

    x_image = nullptr;
    exist_image = nullptr;
    str_magmom = "";

    inputspin.lspin = false;
    inputspin.noncollinear = 0;
    inputspin.time_reversal_symm = 1;
}

void System::deallocate_variables()
{
    if (kdname) {
        deallocate(kdname);
    }
    if (x_image) {
        deallocate(x_image);
    }
    if (is_periodic) {
        deallocate(is_periodic);
    }
    if (exist_image) {
        deallocate(exist_image);
    }
}

void System::set_spin_variables(const size_t nat_in,
                                const bool lspin_in,
                                const int noncol_in,
                                const int trev_sym_in,
                                const double (*magmom_in)[3])
{
    inputspin.lspin = lspin_in;
    inputspin.noncollinear = noncol_in;
    inputspin.time_reversal_symm = trev_sym_in;
    inputspin.magmom.clear();

    std::vector<double> vec(3);
    for (size_t i = 0; i < nat_in; ++i) {
        for (auto j = 0; j < 3; ++j) {
            vec[j] = magmom_in[i][j];
        }
        inputspin.magmom.push_back(vec);
    }
}

const Spin& System::get_spin() const
{
    return inputspin;
}

void System::set_str_magmom(std::string str_magmom_in)
{
    str_magmom = str_magmom_in;
}

const std::string& System::get_str_magmom() const
{
    return str_magmom;
}

const std::vector<std::vector<unsigned int>>& System::get_atomtype_group() const
{
    return atomtype_group;
}

void System::set_atomtype_group()
{
    // In the case of collinear calculation, inputspin moments are considered as scalar
    // variables. Therefore, the same elements with different magnetic moments are
    // considered as different types. In noncollinear calculations,
    // magnetic moments are not considered in this stage. They will be treated
    // separately in symmetry.cpp where inputspin moments will be rotated and flipped
    // using time-reversal symmetry.

    unsigned int i;
    AtomType type_tmp{};
    std::set<AtomType> set_type;
    set_type.clear();

    for (i = 0; i < supercell.number_of_atoms; ++i) {
        type_tmp.element = supercell.kind[i];

        if (spin.noncollinear == 0) {
            type_tmp.magmom = spin.magmom[i][2];
        } else {
            type_tmp.magmom = 0.0;
        }
        set_type.insert(type_tmp);
    }

    const auto natomtypes = set_type.size();
    atomtype_group.resize(natomtypes);

    for (i = 0; i < supercell.number_of_atoms; ++i) {
        int count = 0;
        for (auto it : set_type) {
            if (spin.noncollinear) {
                if (supercell.kind[i] == it.element) {
                    atomtype_group[count].push_back(i);
                }
            } else {
                if ((supercell.kind[i] == it.element)
                    && (std::abs(spin.magmom[i][2] - it.magmom) < eps6)) {
                    atomtype_group[count].push_back(i);
                }
            }
            ++count;
        }
    }
    set_type.clear();
}

void System::generate_coordinate_of_periodic_images()
{
    //
    // Generate Cartesian coordinates of atoms in the neighboring 27 supercells
    //

    unsigned int i;
    int ia, ja, ka;

    const auto nat = supercell.number_of_atoms;
    const auto xf_in = supercell.x_fractional;

    auto icell = 0;
    for (i = 0; i < nat; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            x_image[0][i][j] = xf_in[i][j];
        }
    }
    // Convert to Cartesian coordinate
    frac2cart(x_image[0]);

    for (ia = -1; ia <= 1; ++ia) {
        for (ja = -1; ja <= 1; ++ja) {
            for (ka = -1; ka <= 1; ++ka) {

                if (ia == 0 && ja == 0 && ka == 0) continue;

                ++icell;
                for (i = 0; i < nat; ++i) {
                    x_image[icell][i][0] = xf_in[i][0] + static_cast<double>(ia);
                    x_image[icell][i][1] = xf_in[i][1] + static_cast<double>(ja);
                    x_image[icell][i][2] = xf_in[i][2] + static_cast<double>(ka);
                }
                // Convert to Cartesian coordinate
                frac2cart(x_image[icell]);
            }
        }
    }

    icell = 0;
    exist_image[0] = 1;

    for (ia = -1; ia <= 1; ++ia) {
        for (ja = -1; ja <= 1; ++ja) {
            for (ka = -1; ka <= 1; ++ka) {

                if (ia == 0 && ja == 0 && ka == 0) continue;

                ++icell;

                // When periodic flag is zero along an axis,
                // periodic images along that axis cannot be considered.
                if (((std::abs(ia) == 1) && (is_periodic[0] == 0)) ||
                    ((std::abs(ja) == 1) && (is_periodic[1] == 0)) ||
                    ((std::abs(ka) == 1) && (is_periodic[2] == 0))) {

                    exist_image[icell] = 0;

                } else {

                    exist_image[icell] = 1;
                }
            }
        }
    }
}


void System::print_structure_stdout(const Cell &cell)
{
    using namespace std;
    size_t i;

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
    cout << setw(16) << cell.reciprocal_lattice_vector[0][0];
    cout << setw(15) << cell.reciprocal_lattice_vector[0][1];
    cout << setw(15) << cell.reciprocal_lattice_vector[0][2];
    cout << " : b1" << endl;

    cout << setw(16) << cell.reciprocal_lattice_vector[1][0];
    cout << setw(15) << cell.reciprocal_lattice_vector[1][1];
    cout << setw(15) << cell.reciprocal_lattice_vector[1][2];
    cout << " : b2" << endl;

    cout << setw(16) << cell.reciprocal_lattice_vector[2][0];
    cout << setw(15) << cell.reciprocal_lattice_vector[2][1];
    cout << setw(15) << cell.reciprocal_lattice_vector[2][2];
    cout << " : b3" << endl;
    cout << endl;

    cout << "  Atomic species:" << endl;
    for (i = 0; i < cell.number_of_elems; ++i) {
        cout << setw(6) << i + 1 << setw(5) << kdname[i] << endl;
    }
    cout << endl;

    cout << "  Atomic positions in fractional basis and atomic species" << endl;
    for (i = 0; i < cell.number_of_atoms; ++i) {
        cout << setw(6) << i + 1;
        cout << setw(15) << cell.x_fractional[i][0];
        cout << setw(15) << cell.x_fractional[i][1];
        cout << setw(15) << cell.x_fractional[i][2];
        cout << setw(5) << cell.kind[i] << endl;
    }
    cout << endl << endl;
    cout.unsetf(ios::scientific);
}


void System::print_magmom_stdout() const
{
    using namespace std;

    cout << "  MAGMOM is given. The magnetic moments of each atom are as follows:" << endl;
    for (size_t i = 0; i < inputcell.number_of_atoms; ++i) {
        cout << setw(6) << i + 1;
        cout << setw(5) << inputspin.magmom[i][0];
        cout << setw(5) << inputspin.magmom[i][1];
        cout << setw(5) << inputspin.magmom[i][2];
        cout << endl;
    }
    cout << endl;
    if (inputspin.noncollinear == 0) {
        cout << "  NONCOLLINEAR = 0: magnetic moments are considered as scalar variables." << endl;
    } else if (inputspin.noncollinear == 1) {
        cout << "  NONCOLLINEAR = 1: magnetic moments are considered as vector variables." << endl;
        if (inputspin.time_reversal_symm) {
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
