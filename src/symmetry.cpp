/*
 symmetry.cpp

 Copyright (c) 2014-2018 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "symmetry.h"
#include "error.h"
#include "interaction.h"
#include "memory.h"
#include "mathfunctions.h"
#include "system.h"
#include "timer.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

extern "C" {
#include "spglib.h"
}

using namespace ALM_NS;

Symmetry::Symmetry()
{
    set_default_variables();
}

Symmetry::~Symmetry()
{
    deallocate_variables();
}

void Symmetry::init(const System *system,
                    const int verbosity,
                    Timer *timer)
{
    timer->start_clock("symmetry");

    if (verbosity > 0) {
        std::cout << " SYMMETRY" << std::endl;
        std::cout << " ========" << std::endl << std::endl;
    }


    setup_symmetry_operation(system->get_cell(),
                             system->is_periodic,
                             system->atomtype_group,
                             system->spin,
                             verbosity,
                             SymmData,
                             nsym,
                             nat_prim,
                             ntran,
                             symnum_tran);


    // set_primitive_lattice(system->lavec, system->supercell.number_of_atoms,
    //                       system->kd, system->xcoord,
    //                       lavec_prim, nat_prim,
    //                       kd_prim, xcoord_prim,
    //                       tolerance);

    int nat = system->get_cell().number_of_atoms;

    if (map_sym) {
        deallocate(map_sym);
    }
    allocate(map_sym, nat, nsym);

    if (map_p2s) {
        deallocate(map_p2s);
    }
    allocate(map_p2s, nat_prim, ntran);

    gen_mapping_information(system->get_cell(),
                            system->atomtype_group,
                            SymmData,
                            symnum_tran,
                            map_sym, map_p2s, map_s2p);

    if (verbosity > 0) {
        print_symminfo_stdout();
        timer->print_elapsed();
        std::cout << " -------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
    }


    timer->stop_clock("symmetry");
}

void Symmetry::set_default_variables()
{
    file_sym = "SYMM_INFO";

    // Default values
    nsym = 0;
    printsymmetry = 0;
    map_sym = nullptr;
    map_p2s = nullptr;
    ntran = 0;
    nat_prim = 0;
    tolerance = 1e-6;
    use_internal_symm_finder = false;
}

void Symmetry::deallocate_variables()
{
    if (map_sym) {
        deallocate(map_sym);
    }
    if (map_p2s) {
        deallocate(map_p2s);
    }
}

void Symmetry::setup_symmetry_operation(const Cell &cell,
                                        const int is_periodic[3],
                                        const std::vector<std::vector<unsigned int>> &atomtype_group,
                                        const Spin &spin,
                                        const int verbosity,
                                        std::vector<SymmetryOperation> &SymmData_out,
                                        unsigned int &nsym_out,
                                        unsigned int &nat_prim_out,
                                        unsigned int &ntran_out,
                                        std::vector<int> &symnum_tran_out)
{
    int i, j;

    SymmData_out.clear();

    if (use_internal_symm_finder) {
        findsym_alm(cell, is_periodic, atomtype_group, spin, SymmData_out);
    } else {
        int spgnum;
        std::string spgsymbol;
        findsym_spglib(cell, atomtype_group, spin, tolerance, SymmData_out, spgnum, spgsymbol);

        if (verbosity > 0) {
            std::cout << "  Space group: " << spgsymbol << " (" << std::setw(3) << spgnum << ")" << std::endl;
        }

    }

    // The order in SymmData_out changes for each run because it was generated
    // with OpenMP. Therefore, we sort the list here to have the same result.
    std::sort(SymmData_out.begin() + 1, SymmData_out.end());
    nsym_out = SymmData_out.size();

    if (printsymmetry) {
        std::ofstream ofs_sym;
        if (verbosity > 0) {
            std::cout << "  PRINTSYM = 1: Symmetry information will be stored in SYMM_INFO file."
                << std::endl << std::endl;
        }

        ofs_sym.open(file_sym.c_str(), std::ios::out);
        ofs_sym << nsym << std::endl;

        for (auto &p : SymmData_out) {
            for (i = 0; i < 3; ++i) {
                for (j = 0; j < 3; ++j) {
                    ofs_sym << std::setw(4) << p.rotation[i][j];
                }
            }
            ofs_sym << "  ";
            for (i = 0; i < 3; ++i) {
                ofs_sym << std::setprecision(15) << std::setw(21) << p.tran[i];
            }
            ofs_sym << std::endl;
        }

        ofs_sym.close();
    }

    ntran_out = 0;
    for (i = 0; i < nsym_out; ++i) {
        if (SymmData_out[i].is_translation) ++ntran_out;
    }

    nat_prim_out = cell.number_of_atoms / ntran_out;

    if (cell.number_of_atoms % ntran_out) {
        exit("setup_symmetry_operation",
             "nat != nat_prim * ntran. Something is wrong in the structure.");
    }

    symnum_tran_out.clear();
    for (i = 0; i < nsym_out; ++i) {
        if (SymmData_out[i].is_translation) symnum_tran_out.push_back(i);
    }
}

void Symmetry::findsym_alm(const Cell &cell,
                           const int is_periodic[3],
                           const std::vector<std::vector<unsigned int>> &atomtype_group,
                           const Spin &spin,
                           std::vector<SymmetryOperation> &symop_all)
{
    std::vector<RotationMatrix> LatticeSymmList;

    // Generate rotational matrices that don't change the metric tensor
    LatticeSymmList.clear();
    find_lattice_symmetry(cell.lattice_vector, LatticeSymmList);

    // Generate all the space group operations with translational vectors
    symop_all.clear();
    find_crystal_symmetry(cell,
                          atomtype_group,
                          is_periodic,
                          spin,
                          LatticeSymmList,
                          symop_all);

    LatticeSymmList.clear();
}

void Symmetry::find_lattice_symmetry(const double aa[3][3],
                                     std::vector<RotationMatrix> &LatticeSymmList) const
{
    /*
    Find the rotational matrices that leave the metric tensor invariant.

    Metric tensor G = (g)_{ij} = a_{i} * a_{j} is invariant under crystal symmetry operations T,
    i.e. T^{t}GT = G. Since G can be written as G = A^{t}A, the invariance condition is given by
    (AT)^{t}(AT) = G0 (original).
    */

    int i, j, k;
    int m11, m12, m13, m21, m22, m23, m31, m32, m33;

    int nsym_tmp = 0;
    int mat_tmp[3][3];
    double det, res;
    double rot_tmp[3][3];
    double aa_rot[3][3];

    double metric_tensor[3][3];
    double metric_tensor_rot[3][3];


    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            metric_tensor[i][j] = 0.0;
            for (k = 0; k < 3; ++k) {
                metric_tensor[i][j] += aa[k][i] * aa[k][j];
            }
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            if (i == j) {
                mat_tmp[i][i] = 1;
            } else {
                mat_tmp[i][j] = 0;
            }
        }
    }

    // Identity matrix should be the first entry.
    LatticeSymmList.emplace_back(mat_tmp);

    for (m11 = -1; m11 <= 1; ++m11) {
        for (m12 = -1; m12 <= 1; ++m12) {
            for (m13 = -1; m13 <= 1; ++m13) {
                for (m21 = -1; m21 <= 1; ++m21) {
                    for (m22 = -1; m22 <= 1; ++m22) {
                        for (m23 = -1; m23 <= 1; ++m23) {
                            for (m31 = -1; m31 <= 1; ++m31) {
                                for (m32 = -1; m32 <= 1; ++m32) {
                                    for (m33 = -1; m33 <= 1; ++m33) {

                                        if (m11 == 1 && m12 == 0 && m13 == 0 &&
                                            m21 == 0 && m22 == 1 && m23 == 0 &&
                                            m31 == 0 && m32 == 0 && m33 == 1)
                                            continue;

                                        det = m11 * (m22 * m33 - m32 * m23)
                                            - m21 * (m12 * m33 - m32 * m13)
                                            + m31 * (m12 * m23 - m22 * m13);

                                        if (det != 1 && det != -1) continue;

                                        rot_tmp[0][0] = m11;
                                        rot_tmp[0][1] = m12;
                                        rot_tmp[0][2] = m13;
                                        rot_tmp[1][0] = m21;
                                        rot_tmp[1][1] = m22;
                                        rot_tmp[1][2] = m23;
                                        rot_tmp[2][0] = m31;
                                        rot_tmp[2][1] = m32;
                                        rot_tmp[2][2] = m33;

                                        // Here, aa_rot = aa * rot_tmp is correct.
                                        matmul3(aa_rot, aa, rot_tmp);

                                        for (i = 0; i < 3; ++i) {
                                            for (j = 0; j < 3; ++j) {
                                                metric_tensor_rot[i][j] = 0.0;
                                                for (k = 0; k < 3; ++k) {
                                                    metric_tensor_rot[i][j] += aa_rot[k][i] * aa_rot[k][j];
                                                }
                                            }
                                        }

                                        res = 0.0;
                                        for (i = 0; i < 3; ++i) {
                                            for (j = 0; j < 3; ++j) {
                                                res += std::pow(metric_tensor[i][j] - metric_tensor_rot[i][j], 2.0);
                                            }
                                        }

                                        // Metric tensor is invariant under symmetry operations.

                                        if (res < tolerance * tolerance) {
                                            ++nsym_tmp;
                                            for (i = 0; i < 3; ++i) {
                                                for (j = 0; j < 3; ++j) {
                                                    mat_tmp[i][j] = static_cast<int>(rot_tmp[i][j]);
                                                }
                                            }
                                            LatticeSymmList.emplace_back(mat_tmp);
                                        }

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (LatticeSymmList.size() > 48) {
        exit("find_lattice_symmetry", "Number of lattice symmetry is larger than 48.");
    }
}

void Symmetry::find_crystal_symmetry(const Cell &cell,
                                     const std::vector<std::vector<unsigned int>> &atomtype_group,
                                     const int is_periodic[3],
                                     const Spin &spin,
                                     const std::vector<RotationMatrix> &LatticeSymmList,
                                     std::vector<SymmetryOperation> &CrystalSymmList)
{
    unsigned int i, j;
    unsigned int iat, jat, kat, lat;
    double x_rot[3], x_tmp[3];
    double rot[3][3], rot_tmp[3][3], rot_cart[3][3];
    double mag[3], mag_rot[3];
    double tran[3];
    double x_rot_tmp[3];
    double tmp[3];
    double diff;
    auto nclass = atomtype_group.size();

    int rot_int[3][3];

    int ii, jj, kk;
    unsigned int itype;

    bool is_found;
    bool isok;
    bool mag_sym1, mag_sym2;
    bool is_identity_matrix;


    // Add identity matrix first.
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            if (i == j) {
                rot_int[i][j] = 1;
                rot_cart[i][j] = 1.0;
            } else {
                rot_int[i][j] = 0;
                rot_cart[i][j] = 0.0;
            }
        }
        tran[i] = 0.0;
    }

    CrystalSymmList.emplace_back(rot_int,
                                 tran,
                                 rot_cart,
                                 is_compatible(rot_int),
                                 is_compatible(rot_cart),
                                 is_translation(rot_int));

    for (auto &it_latsym : LatticeSymmList) {

        iat = atomtype_group[0][0];

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                rot[i][j] = static_cast<double>(it_latsym.mat[i][j]);
            }
        }

        for (i = 0; i < 3; ++i) x_tmp[i] = cell.x_fractional[iat][i];
        rotvec(x_rot, x_tmp, rot);

#ifdef _OPENMP
#pragma omp parallel for private(jat, tran, isok, kat, x_tmp, x_rot_tmp, is_found, lat, tmp, diff, \
    i, j, itype, jj, kk, is_identity_matrix, mag, mag_rot, rot_tmp, rot_cart, mag_sym1, mag_sym2)
#endif
        for (ii = 0; ii < atomtype_group[0].size(); ++ii) {
            jat = atomtype_group[0][ii];

            for (i = 0; i < 3; ++i) {
                tran[i] = cell.x_fractional[jat][i] - x_rot[i];
                tran[i] = tran[i] - nint(tran[i]);
            }

            if ((std::abs(tran[0]) > eps12 && !is_periodic[0]) ||
                (std::abs(tran[1]) > eps12 && !is_periodic[1]) ||
                (std::abs(tran[2]) > eps12 && !is_periodic[2]))
                continue;

            is_identity_matrix =
            (std::pow(rot[0][0] - 1.0, 2) + std::pow(rot[0][1], 2) + std::pow(rot[0][2], 2)
                + std::pow(rot[1][0], 2) + std::pow(rot[1][1] - 1.0, 2) + std::pow(rot[1][2], 2)
                + std::pow(rot[2][0], 2) + std::pow(rot[2][1], 2) + std::pow(rot[2][2] - 1.0, 2)
                + std::pow(tran[0], 2) + std::pow(tran[1], 2) + std::pow(tran[2], 2)) < eps12;
            if (is_identity_matrix) continue;

            isok = true;

            for (itype = 0; itype < nclass; ++itype) {

                for (jj = 0; jj < atomtype_group[itype].size(); ++jj) {

                    kat = atomtype_group[itype][jj];

                    for (i = 0; i < 3; ++i) x_tmp[i] = cell.x_fractional[kat][i];
                    rotvec(x_rot_tmp, x_tmp, rot);

                    for (i = 0; i < 3; ++i) {
                        x_rot_tmp[i] += tran[i];
                    }

                    is_found = false;

                    for (kk = 0; kk < atomtype_group[itype].size(); ++kk) {

                        lat = atomtype_group[itype][kk];

                        for (i = 0; i < 3; ++i) {
                            tmp[i] = std::fmod(std::abs(cell.x_fractional[lat][i] - x_rot_tmp[i]), 1.0);
                            tmp[i] = std::min<double>(tmp[i], 1.0 - tmp[i]);
                        }
                        diff = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
                        if (diff < tolerance * tolerance) {
                            is_found = true;
                            break;
                        }
                    }

                    if (!is_found) isok = false;
                }
            }

            if (isok) {
                matmul3(rot_tmp, rot, cell.reciprocal_lattice_vector);
                matmul3(rot_cart, cell.lattice_vector, rot_tmp);

                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        rot_cart[i][j] /= (2.0 * pi);
                    }
                }

                if (spin.lspin && spin.noncollinear) {

                    for (i = 0; i < 3; ++i) {
                        mag[i] = spin.magmom[jat][i];
                        mag_rot[i] = spin.magmom[iat][i];
                    }

                    rotvec(mag_rot, mag_rot, rot_cart);

                    // In the case of improper rotation, the factor -1 should be multiplied
                    // because the inversion operation doesn't flip the spin.
                    if (!is_proper(rot_cart)) {
                        for (i = 0; i < 3; ++i) {
                            mag_rot[i] = -mag_rot[i];
                        }
                    }

                    mag_sym1 = (std::pow(mag[0] - mag_rot[0], 2.0)
                        + std::pow(mag[1] - mag_rot[1], 2.0)
                        + std::pow(mag[2] - mag_rot[2], 2.0)) < eps6;

                    mag_sym2 = (std::pow(mag[0] + mag_rot[0], 2.0)
                        + std::pow(mag[1] + mag_rot[1], 2.0)
                        + std::pow(mag[2] + mag_rot[2], 2.0)) < eps6;

                    if (!mag_sym1 && !mag_sym2) {
                        isok = false;
                    } else if (!mag_sym1 && mag_sym2 && !spin.time_reversal_symm) {
                        isok = false;
                    }
                }
            }


            if (isok) {
#ifdef _OPENMP
#pragma omp critical
#endif
                CrystalSymmList.emplace_back(it_latsym.mat,
                                             tran,
                                             rot_cart,
                                             is_compatible(it_latsym.mat),
                                             is_compatible(rot_cart),
                                             is_translation(it_latsym.mat));
            }
        }

    }
}

void Symmetry::findsym_spglib(const Cell &cell,
                              const std::vector<std::vector<unsigned int>> &atomtype_group,
                              const Spin &spin,
                              const double symprec,
                              std::vector<SymmetryOperation> &symop_all,
                              int &spgnum,
                              std::string &spgsymbol)
{
    int i, j;
    double (*position)[3];
    double (*translation)[3];
    int (*rotation)[3][3];
    char symbol[11];
    double aa_tmp[3][3];
    int *types_tmp;

    auto nat = cell.number_of_atoms;

    if (spin.lspin && spin.noncollinear) {
        exit("findsym_spglib", "NONCOLLINEAR spin is not supported in spglib.");
    }

    allocate(position, nat);
    allocate(types_tmp, nat);


    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            aa_tmp[i][j] = cell.lattice_vector[i][j];
        }
    }
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            position[i][j] = cell.x_fractional[i][j];
        }
    }

    if (spin.lspin) {
        for (i = 0; i < atomtype_group.size(); ++i) {
            for (j = 0; j < atomtype_group[i].size(); ++j) {
                types_tmp[atomtype_group[i][j]] = i;
            }
        }
    } else {
        for (i = 0; i < nat; ++i) {
            types_tmp[i] = cell.kind[i];
        }
    }

    // First find the number of symmetry operations
    auto nsym = spg_get_multiplicity(aa_tmp, position, types_tmp, nat, symprec);

    if (nsym == 0) exit("findsym_spglib", "Error occured in spg_get_multiplicity");

    allocate(translation, nsym);
    allocate(rotation, nsym);

    // Store symmetry operations
    nsym = spg_get_symmetry(rotation, translation, nsym,
                            aa_tmp, position, types_tmp, nat, symprec);

    spgnum = spg_get_international(symbol, aa_tmp, position, types_tmp, nat, symprec);
    spgsymbol = std::string(symbol);

    // Copy symmetry information
    symop_all.clear();
    double rot_cartesian[3][3];

    for (i = 0; i < nsym; ++i) {

        symop_in_cart(rot_cartesian,
                      rotation[i],
                      cell.lattice_vector,
                      cell.reciprocal_lattice_vector);

        symop_all.emplace_back(rotation[i],
                               translation[i],
                               rot_cartesian,
                               is_compatible(rotation[i]),
                               is_compatible(rot_cartesian),
                               is_translation(rotation[i]));

    }


    deallocate(rotation);
    deallocate(translation);
    deallocate(types_tmp);
    deallocate(position);
}

void Symmetry::symop_in_cart(double rot_cart[3][3],
                             const int rot_lattice[3][3],
                             const double lavec[3][3],
                             const double rlavec[3][3]) const
{
    int i, j;
    double sym_tmp[3][3];
    double tmp[3][3];

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            sym_tmp[i][j] = static_cast<double>(rot_lattice[i][j]);
        }
    }

    matmul3(tmp, sym_tmp, rlavec);
    matmul3(rot_cart, lavec, tmp);

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            rot_cart[i][j] = rot_cart[i][j] / (2.0 * pi);
        }
    }
}


void Symmetry::print_symminfo_stdout() const
{
    int i;

    std::cout << "  Number of symmetry operations = " << SymmData.size() << std::endl;
    std::cout << std::endl;
    if (ntran > 1) {
        std::cout << "  Given system is not a primitive cell." << std::endl;
        std::cout << "  There are " << std::setw(5)
            << ntran << " translation operations." << std::endl;
    } else {
        std::cout << "  Given system is a primitive cell." << std::endl;
    }
    std::cout << "  Primitive cell contains " << nat_prim << " atoms" << std::endl;

    std::cout << std::endl;
    std::cout << "  **Cell-Atom Correspondens Below**" << std::endl;
    std::cout << std::setw(6) << " CELL" << " | " << std::setw(5) << "ATOM" << std::endl;

    for (int i = 0; i < ntran; ++i) {
        std::cout << std::setw(6) << i + 1 << " | ";
        for (int j = 0; j < nat_prim; ++j) {
            std::cout << std::setw(5) << map_p2s[j][i] + 1;
            if ((j + 1) % 5 == 0) {
                std::cout << std::endl << "       | ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Symmetry::gen_mapping_information(const Cell &cell,
                                       const std::vector<std::vector<unsigned int>> &atomtype_group,
                                       const std::vector<SymmetryOperation> &symmdata,
                                       const std::vector<int> &symnum_tran,
                                       int **map_sym,
                                       int **map_p2s,
                                       std::vector<Maps> &map_s2p) const
{
    int isym, iat, jat;
    int i, j;
    int itype;
    int ii, jj;
    double xnew[3], x_tmp[3];
    double tmp[3], diff;
    double rot_double[3][3];

    for (iat = 0; iat < cell.number_of_atoms; ++iat) {
        for (isym = 0; isym < symmdata.size(); ++isym) {
            map_sym[iat][isym] = -1;
        }
    }
    int nsym = symmdata.size();

    // This part may be incompatible with the tolerance used in spglib
    auto natomtypes = atomtype_group.size();

#ifdef _OPENMP
#pragma omp parallel for private(i, j, rot_double, itype, ii, iat, x_tmp, xnew, jj, jat, tmp, diff, isym)
#endif
    for (isym = 0; isym < nsym; ++isym) {

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                rot_double[i][j] = static_cast<double>(symmdata[isym].rotation[i][j]);
            }
        }

        for (itype = 0; itype < natomtypes; ++itype) {

            for (ii = 0; ii < atomtype_group[itype].size(); ++ii) {

                iat = atomtype_group[itype][ii];

                for (i = 0; i < 3; ++i) x_tmp[i] = cell.x_fractional[iat][i];
                rotvec(xnew, x_tmp, rot_double);

                for (i = 0; i < 3; ++i) xnew[i] += symmdata[isym].tran[i];

                for (jj = 0; jj < atomtype_group[itype].size(); ++jj) {

                    jat = atomtype_group[itype][jj];

                    for (i = 0; i < 3; ++i) {
                        tmp[i] = std::fmod(std::abs(cell.x_fractional[jat][i] - xnew[i]), 1.0);
                        tmp[i] = std::min<double>(tmp[i], 1.0 - tmp[i]);
                    }
                    diff = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
                    if (diff < tolerance * tolerance) {
                        map_sym[iat][isym] = jat;
                        break;
                    }
                }
                if (map_sym[iat][isym] == -1) {
                    exit("gen_mapping_information",
                         "cannot find symmetry for operation # ",
                         isym + 1);
                }
            }
        }
    }

    double shift[3];
    double pos[3], pos_std[3];

    // Generate map_p2s (primitive --> super)

    bool *is_checked;
    allocate(is_checked, cell.number_of_atoms);

    for (i = 0; i < cell.number_of_atoms; ++i) is_checked[i] = false;

    jat = 0;
    int ntran = symnum_tran.size();
    int atomnum_translated;
    for (iat = 0; iat < cell.number_of_atoms; ++iat) {

        if (is_checked[iat]) continue;
        for (i = 0; i < ntran; ++i) {
            atomnum_translated = map_sym[iat][symnum_tran[i]];
            map_p2s[jat][i] = atomnum_translated;
            is_checked[atomnum_translated] = true;
        }
        ++jat;
    }

    deallocate(is_checked);

    // Generate map_s2p (super --> primitive)

    int nat_prim = cell.number_of_atoms / ntran;
    map_s2p.clear();
    map_s2p.resize(cell.number_of_atoms);

    for (iat = 0; iat < nat_prim; ++iat) {
        for (i = 0; i < ntran; ++i) {
            atomnum_translated = map_p2s[iat][i];
            map_s2p[atomnum_translated].atom_num = iat;
            map_s2p[atomnum_translated].tran_num = i;
        }
    }
}

bool Symmetry::is_translation(const int rot[3][3]) const
{
    bool ret =
        rot[0][0] == 1 && rot[0][1] == 0 && rot[0][2] == 0 &&
        rot[1][0] == 0 && rot[1][1] == 1 && rot[1][2] == 0 &&
        rot[2][0] == 0 && rot[2][1] == 0 && rot[2][2] == 1;

    return ret;
}

template <typename T>
bool Symmetry::is_compatible(const T rot[3][3],
                             const double tolerance_zero)
{
    int nfinite = 0;
    double rot_double[3][3];

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            rot_double[i][j] = static_cast<double>(rot[i][j]);
            if (std::abs(rot_double[i][j]) > tolerance_zero) ++nfinite;
        }
    }

    return (nfinite == 3);
}


bool Symmetry::is_proper(const double rot[3][3]) const
{
    double det = rot[0][0] * (rot[1][1] * rot[2][2] - rot[2][1] * rot[1][2])
        - rot[1][0] * (rot[0][1] * rot[2][2] - rot[2][1] * rot[0][2])
        + rot[2][0] * (rot[0][1] * rot[1][2] - rot[1][1] * rot[0][2]);

    if (std::abs(det - 1.0) < eps12) {
        return true;
    }
    if (std::abs(det + 1.0) < eps12) {
        return false;
    }
    exit("is_proper", "This cannot happen.");
}

void Symmetry::set_primitive_lattice(const double aa[3][3],
                                     const int nat,
                                     const int *kd,
                                     double **x,
                                     double aa_prim[3][3],
                                     unsigned int &nat_prim,
                                     int *kd_prim,
                                     double **x_prim,
                                     const double symprec) const
{
    int i, j;
    int *types_tmp;
    double (*position)[3];

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            aa_prim[i][j] = aa[i][j];
        }
    }

    allocate(position, nat);
    allocate(types_tmp, nat);

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            position[i][j] = x[i][j];
        }
        types_tmp[i] = kd[i];
    }

    //    nat_prim = spg_find_primitive(aa_prim, position, types_tmp, nat, symprec);
    nat_prim = spg_standardize_cell(aa_prim,
                                    position,
                                    types_tmp,
                                    nat, 1, 0,
                                    symprec);

    for (i = 0; i < nat_prim; ++i) {
        for (j = 0; j < 3; ++j) {
            x_prim[i][j] = position[i][j];
        }
        kd_prim[i] = types_tmp[i];
    }

    deallocate(position);
    deallocate(types_tmp);
}
