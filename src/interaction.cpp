/*
interaction.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "interaction.h"
#include "combination.h"
#include "constants.h"
#include "fcs.h"
#include "mathfunctions.h"
#include "memory.h"
#include "symmetry.h"
#include "system.h"
#include "timer.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>

using namespace ALM_NS;

Interaction::Interaction()
{
    set_default_variables();
}

Interaction::~Interaction()
{
    deallocate_variables();
}

void Interaction::init(ALM *alm)
{
    alm->timer->start_clock("interaction");

    int i, j, k;
    int nat = alm->system->supercell.number_of_atoms;
    int nkd = alm->system->supercell.number_of_elems;


    std::cout << " INTERACTION" << std::endl;
    std::cout << " ===========" << std::endl << std::endl;

    set_ordername();

    if (distall) {
        deallocate(distall);
    }
    allocate(distall, nat, nat);

    if (mindist_pairs) {
        deallocate(mindist_pairs);
    }
    allocate(mindist_pairs, nat, nat);

    if (interaction_pair) {
        deallocate(interaction_pair);
    }
    allocate(interaction_pair, maxorder, alm->symmetry->nat_prim);

    if (interaction_cluster) {
        deallocate(interaction_cluster);
    }
    allocate(interaction_cluster, maxorder, alm->symmetry->nat_prim);

    if (cluster_list) {
        deallocate(cluster_list);
    }
    allocate(cluster_list, maxorder);

    get_pairs_of_minimum_distance(nat,
                                  alm->system->x_image,
                                  alm->system->exist_image);

    set_interaction_by_cutoff(alm->system->supercell.number_of_atoms,
                              alm->system->supercell.kind,
                              alm->symmetry->nat_prim,
                              alm->symmetry->map_p2s,
                              cutoff_radii,
                              interaction_pair);

    calc_interaction_clusters(alm->symmetry->nat_prim,
                              alm->system->supercell.kind,
                              alm->symmetry->map_p2s,
                              interaction_pair,
                              alm->system->x_image,
                              alm->system->exist_image,
                              interaction_cluster);

    generate_pairs(alm->symmetry->nat_prim,
                   alm->symmetry->map_p2s,
                   cluster_list,
                   interaction_cluster);


    std::cout << "  +++ Cutoff Radii Matrix in Bohr Unit (NKD x NKD matrix) +++" << std::endl;

    for (i = 0; i < maxorder; ++i) {
        std::cout << "  " << std::setw(9) << str_order[i] << std::endl;
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                if (cutoff_radii[i][j][k] < 0.0) {
                    std::cout << std::setw(9) << "None";
                } else {
                    std::cout << std::setw(9) << cutoff_radii[i][j][k];
                }
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    print_interaction_information(alm->symmetry->nat_prim,
                                  alm->symmetry->map_p2s,
                                  alm->system->supercell.kind,
                                  alm->system->kdname,
                                  interaction_pair);

    print_neighborlist(alm->system->supercell.number_of_atoms,
                       alm->symmetry->nat_prim,
                       alm->symmetry->map_p2s,
                       alm->system->supercell.kind,
                       alm->system->kdname);

    for (i = 0; i < maxorder; ++i) {
        if (i + 2 > nbody_include[i]) {
            std::cout << "  For " << std::setw(8) << str_order[i] << ", ";
            std::cout << "interactions related to more than" << std::setw(2) << nbody_include[i];
            std::cout << " atoms will be neglected." << std::endl;
        }
    }


    alm->timer->print_elapsed();
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    alm->timer->stop_clock("interaction");
}

void Interaction::generate_pairs(const int natmin,
                                 int **map_p2s,
                                 std::set<IntList> *pair_out,
                                 std::set<InteractionCluster> **interaction_cluster)
{
    int i, j;
    int iat;
    int *pair_tmp;

    for (int order = 0; order < maxorder; ++order) {

        pair_out[order].clear();

        allocate(pair_tmp, order + 2);

        for (i = 0; i < natmin; ++i) {

            iat = map_p2s[i][0];

            for (auto it = interaction_cluster[order][i].begin();
                 it != interaction_cluster[order][i].end(); ++it) {

                pair_tmp[0] = iat;
                for (j = 0; j < order + 1; ++j) {
                    pair_tmp[j + 1] = (*it).atom[j];
                }

                insort(order + 2, pair_tmp);

                // Ignore many-body case 
                // if (!satisfy_nbody_rule(order + 2, pair_tmp, order)) continue;
                pair_out[order].insert(IntList(order + 2, pair_tmp));
            }
        }
        deallocate(pair_tmp);
    }
}

void Interaction::set_default_variables()
{
    maxorder = 0;
    nbody_include = nullptr;
    cutoff_radii = nullptr;
    distall = nullptr;
    mindist_pairs = nullptr;
    cluster_list = nullptr;
    interaction_pair = nullptr;
    interaction_cluster = nullptr;
}

void Interaction::deallocate_variables()
{
    if (cluster_list) {
        deallocate(cluster_list);
    }
    if (mindist_pairs) {
        deallocate(mindist_pairs);
    }
    if (interaction_pair) {
        deallocate(interaction_pair);
    }
    if (interaction_cluster) {
        deallocate(interaction_cluster);
    }
    if (distall) {
        deallocate(distall);
    }
}

double Interaction::distance(double *x1, double *x2)
{
    double dist;
    dist = std::pow(x1[0] - x2[0], 2) + std::pow(x1[1] - x2[1], 2) + std::pow(x1[2] - x2[2], 2);
    dist = std::sqrt(dist);

    return dist;
}

void Interaction::get_pairs_of_minimum_distance(int nat,
                                                double ***xc_in,
                                                int *exist)
{
    //
    // Calculate the minimum distance between atom i and j 
    // under the periodic boundary conditions
    //
    int icell;
    int i, j, k;
    double dist_tmp;
    double vec[3];


    for (i = 0; i < nat; ++i) {
        for (j = 0; j < nat; ++j) {

            distall[i][j].clear();

            for (icell = 0; icell < 27; ++icell) {

                if (exist[icell]) {

                    dist_tmp = distance(xc_in[0][i], xc_in[icell][j]);

                    for (k = 0; k < 3; ++k) vec[k] = xc_in[icell][j][k] - xc_in[0][i][k];

                    distall[i][j].emplace_back(DistInfo(icell, dist_tmp, vec));
                }
            }
            std::sort(distall[i][j].begin(), distall[i][j].end());
        }
    }

    // Construct pairs of minimum distance.

    double dist_min;
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < nat; ++j) {
            mindist_pairs[i][j].clear();

            dist_min = distall[i][j][0].dist;
            for (auto it = distall[i][j].cbegin(); it != distall[i][j].cend(); ++it) {
                // The tolerance below (1.e-3) should be chosen so that 
                // the mirror images with equal distances are found correctly.
                // If this fails, the phonon dispersion would be incorrect.
                if (std::abs((*it).dist - dist_min) < 1.0e-3) {
                    mindist_pairs[i][j].emplace_back(DistInfo(*it));
                }
            }
        }
    }
}

void Interaction::print_neighborlist(const int nat,
                                     const int natmin,
                                     int **map_p2s,
                                     const std::vector<int> &kd,
                                     std::string *kdname)
{
    //
    // Print the list of neighboring atoms and distances
    //
    int i, j, k;
    int iat;
    int icount;

    double dist_tmp;
    std::vector<DistList> *neighborlist;

    allocate(neighborlist, natmin);

    for (i = 0; i < natmin; ++i) {
        neighborlist[i].clear();

        iat = map_p2s[i][0];

        for (j = 0; j < nat; ++j) {
            neighborlist[i].emplace_back(DistList(j, mindist_pairs[iat][j][0].dist));
        }
        std::sort(neighborlist[i].begin(), neighborlist[i].end());
    }

    std::cout << std::endl;
    std::cout << "  List of neighboring atoms below." << std::endl;
    std::cout << "  Format [N th-nearest shell, distance in Bohr (Number of atoms on the shell)]"
        << std::endl << std::endl;

    int nthnearest;
    std::vector<int> atomlist;

    for (i = 0; i < natmin; ++i) {

        nthnearest = 0;
        atomlist.clear();

        iat = map_p2s[i][0];
        std::cout << std::setw(5) << iat + 1 << " ("
            << std::setw(3) << kdname[kd[iat] - 1] << "): ";

        dist_tmp = 0.0;

        for (j = 0; j < nat; ++j) {

            if (neighborlist[i][j].dist < eps8) continue; // distance is zero

            if (std::abs(neighborlist[i][j].dist - dist_tmp) > eps6) {

                if (!atomlist.empty()) {
                    nthnearest += 1;

                    if (nthnearest > 1) std::cout << std::setw(13) << " ";

                    std::cout << std::setw(3) << nthnearest << std::setw(10) << dist_tmp
                        << " (" << std::setw(3) << atomlist.size() << ") -";

                    icount = 0;
                    for (k = 0; k < atomlist.size(); ++k) {

                        if (icount % 4 == 0 && icount > 0) {
                            std::cout << std::endl;
                            std::cout << std::setw(34) << " ";
                        }
                        ++icount;

                        std::cout << std::setw(4) << atomlist[k] + 1;
                        std::cout << "(" << std::setw(3)
                            << kdname[kd[atomlist[k]] - 1] << ")";

                    }
                    std::cout << std::endl;
                }


                dist_tmp = neighborlist[i][j].dist;
                atomlist.clear();
                atomlist.push_back(neighborlist[i][j].atom);
            } else {
                atomlist.push_back(neighborlist[i][j].atom);
            }
        }
        if (!atomlist.empty()) {
            nthnearest += 1;

            if (nthnearest > 1) std::cout << std::setw(13) << " ";

            std::cout << std::setw(3) << nthnearest << std::setw(10) << dist_tmp
                << " (" << std::setw(3) << atomlist.size() << ") -";

            icount = 0;
            for (k = 0; k < atomlist.size(); ++k) {

                if (icount % 4 == 0 && icount > 0) {
                    std::cout << std::endl;
                    std::cout << std::setw(34) << " ";
                }
                ++icount;

                std::cout << std::setw(4) << atomlist[k] + 1;
                std::cout << "(" << std::setw(3)
                    << kdname[kd[atomlist[k]] - 1] << ")";

            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    atomlist.clear();
    deallocate(neighborlist);
}


void Interaction::generate_interaction_information_by_cutoff(const int nat,
                                                             const int natmin,
                                                             const std::vector<int> &kd,
                                                             int **map_p2s,
                                                             double **rc,
                                                             std::vector<int> *interaction_list)
{
    int i, iat, jat, ikd, jkd;
    double cutoff_tmp;

    for (i = 0; i < natmin; ++i) {

        interaction_list[i].clear();

        iat = map_p2s[i][0];
        ikd = kd[iat] - 1;

        for (jat = 0; jat < nat; ++jat) {

            jkd = kd[jat] - 1;

            cutoff_tmp = rc[ikd][jkd];

            if (cutoff_tmp < 0.0) {

                // Cutoff 'None'
                interaction_list[i].push_back(jat);

            } else {

                if (mindist_pairs[iat][jat][0].dist <= cutoff_tmp) {
                    interaction_list[i].push_back(jat);
                }
            }
        }
    }
}

void Interaction::set_interaction_by_cutoff(const unsigned int nat,
                                            const std::vector<int> &kd,
                                            const unsigned int nat_prim,
                                            int **map_p2s_in,
                                            double ***rcs,
                                            std::vector<int> **interaction_pair_out)
{
    for (int order = 0; order < maxorder; ++order) {
        generate_interaction_information_by_cutoff(nat,
                                                   nat_prim,
                                                   kd,
                                                   map_p2s_in,
                                                   rcs[order],
                                                   interaction_pair_out[order]);
    }
}

void Interaction::print_interaction_information(const int natmin,
                                                int **map_p2s,
                                                const std::vector<int> &kd,
                                                std::string *kdname,
                                                std::vector<int> **interaction_list)
{
    int order;
    int i, iat;
    std::vector<int> intlist;

    std::cout << std::endl;
    std::cout << "  List of interacting atom pairs considered for each order:" << std::endl;

    for (order = 0; order < maxorder; ++order) {

        std::cout << std::endl << "   ***" << str_order[order] << "***" << std::endl;

        for (i = 0; i < natmin; ++i) {

            if (interaction_list[order][i].empty()) {
                std::cout << "   No interacting atoms! Skipped." << std::endl;
                continue; // no interaction
            }

            iat = map_p2s[i][0];

            intlist.clear();
            for (auto &it : interaction_list[order][i]) {
                intlist.push_back(it);
            }
            std::sort(intlist.begin(), intlist.end());

            // write atoms inside the cutoff radius
            int id = 0;
            std::cout << "    Atom " << std::setw(5) << iat + 1
                << "(" << std::setw(3) << kdname[kd[iat] - 1]
                << ")" << " interacts with atoms ... " << std::endl;

            for (int id = 0; id < intlist.size(); ++id) {
                if (id % 6 == 0) {
                    if (id == 0) {
                        std::cout << "   ";
                    } else {
                        std::cout << std::endl;
                        std::cout << "   ";
                    }
                }
                std::cout << std::setw(5) << intlist[id] + 1 << "("
                    << std::setw(3) << kdname[kd[intlist[id]] - 1] << ")";
            }

            std::cout << std::endl << std::endl;
            std::cout << "    Number of total interaction pairs = "
                << interaction_list[order][i].size() << std::endl << std::endl;
        }
    }

    std::cout << std::endl;
}


bool Interaction::is_incutoff(const int n,
                              int *atomnumlist,
                              const int order,
                              const std::vector<int> &kd)
{
    int iat, jat;
    int ikd, jkd;
    double cutoff_tmp;

    for (int i = 0; i < n; ++i) {
        iat = atomnumlist[i];
        ikd = kd[iat] - 1;

        for (int j = i + 1; j < n; ++j) {
            jat = atomnumlist[j];
            jkd = kd[jat] - 1;

            cutoff_tmp = cutoff_radii[order][ikd][jkd];

            if (cutoff_tmp >= 0.0 &&
                (mindist_pairs[iat][jat][0].dist > cutoff_tmp)) {
                return false;
            }

        }
    }
    return true;
}


void Interaction::set_ordername()
{
    str_order.resize(maxorder);
    str_order[0] = "HARMONIC";

    for (int i = 1; i < maxorder; ++i) {
        str_order[i] = "ANHARM" + std::to_string(i + 2);
    }
}

int Interaction::nbody(const int n, const int *arr)
{
    std::vector<int> v(n);

    for (unsigned int i = 0; i < n; ++i) {
        v[i] = arr[i];
    }
    std::stable_sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());

    int ret = v.size();
    v.clear();

    return ret;
}

bool Interaction::satisfy_nbody_rule(const int nelem,
                                     const int *arr,
                                     const int order)
{
    return nbody(nelem, arr) <= nbody_include[order];
}


void Interaction::calc_interaction_clusters(const int natmin,
                                            const std::vector<int> &kd,
                                            int **map_p2s,
                                            std::vector<int> **interaction_pair_in,
                                            double ***x_image,
                                            int *exist,
                                            std::set<InteractionCluster> **interaction_cluster_out)
{
    //
    // Calculate the complete set of interaction clusters for all orders.
    //

    for (int order = 0; order < maxorder; ++order) {
        set_interaction_cluster(order,
                                natmin,
                                kd,
                                map_p2s,
                                interaction_pair_in[order],
                                x_image,
                                exist,
                                interaction_cluster_out[order]);

    }
}


void Interaction::set_interaction_cluster(const int order,
                                          const int natmin,
                                          const std::vector<int> &kd,
                                          int **map_p2s,
                                          std::vector<int> *interaction_pair_in,
                                          double ***x_image,
                                          int *exist,
                                          std::set<InteractionCluster> *interaction_cluster_out)
{
    //
    // Calculate a set of interaction clusters for the given order
    //

    int j, k;
    int iat, jat;
    int ii;
    int ikd, jkd;
    int icount;
    int idata;
    unsigned int ielem;

    double dist_tmp, rc_tmp;
    double distmax;

    bool isok;

    int *list_now;

    std::vector<int> cell_vector;
    std::vector<double> dist_vector;
    std::vector<std::vector<int>> pairs_icell, comb_cell, comb_cell_min;
    std::vector<std::vector<int>> comb_cell_atom_center;
    std::vector<int> accum_tmp;
    std::vector<int> atom_tmp, cell_tmp;
    std::vector<int> intpair_uniq, cellpair;
    std::vector<int> group_atom;
    std::vector<int> data_now;
    std::vector<MinDistList> distance_list;
    std::vector<std::vector<int>> data_vec;

    allocate(list_now, order + 2);

    for (int i = 0; i < natmin; ++i) {

        interaction_cluster_out[i].clear();

        iat = map_p2s[i][0];
        ikd = kd[iat] - 1;
        list_now[0] = iat;

        // List of 2-body interaction pairs
        std::vector<int> intlist(interaction_pair_in[i]);
        std::sort(intlist.begin(), intlist.end()); // Need to sort here

        if (order == 0) {

            // Harmonic term

            for (ielem = 0; ielem < intlist.size(); ++ielem) {

                jat = intlist[ielem];
                list_now[1] = jat;

                if (!satisfy_nbody_rule(2, list_now, 0)) continue;

                comb_cell_min.clear();
                atom_tmp.clear();
                atom_tmp.push_back(jat);

                for (j = 0; j < mindist_pairs[iat][jat].size(); ++j) {
                    cell_tmp.clear();
                    cell_tmp.push_back(mindist_pairs[iat][jat][j].cell);
                    comb_cell_min.push_back(cell_tmp);
                }
                distmax = mindist_pairs[iat][jat][0].dist;
                interaction_cluster_out[i].insert(InteractionCluster(atom_tmp,
                                                                     comb_cell_min,
                                                                     distmax));
            }

        } else if (order > 0) {

            // Anharmonic terms

            data_vec.clear();

            // First, we generate all possible combinations of interaction clusters.
            CombinationWithRepetition<int> g(intlist.begin(), intlist.end(), order + 1);
            do {
                std::vector<int> data = g.now();

                list_now[0] = iat;
                for (j = 0; j < order + 1; ++j) list_now[j + 1] = data[j];

                // Save as a candidate if the cluster satisfies the NBODY-rule.
                if (satisfy_nbody_rule(order + 2, list_now, order)) {
                    data_vec.emplace_back(data);
                }

            } while (g.next());

            intlist.clear();

            int ndata = data_vec.size();

            for (idata = 0; idata < ndata; ++idata) {

                data_now = data_vec[idata];

                // Uniq the list of atoms in data like as follows:
                // cubic   term : (i, i) --> (i) x 2
                // quartic term : (i, i, j) --> (i, j) x (2, 1)
                intpair_uniq.clear();
                group_atom.clear();
                icount = 1;

                for (j = 0; j < order; ++j) {
                    if (data_now[j] == data_now[j + 1]) {
                        ++icount;
                    } else {
                        group_atom.push_back(icount);
                        intpair_uniq.push_back(data_now[j]);
                        icount = 1;
                    }
                }
                group_atom.push_back(icount);
                intpair_uniq.push_back(data_now[order]);

                pairs_icell.clear();
                for (j = 0; j < intpair_uniq.size(); ++j) {
                    jat = intpair_uniq[j];
                    jkd = kd[jat] - 1;

                    rc_tmp = cutoff_radii[order][ikd][jkd];
                    cell_vector.clear();

                    // Loop over the cell images of atom 'jat' and add to the list 
                    // as a candidate for the interaction cluster.
                    // The mirror images whose distance is larger than the minimum value
                    // of the distance(iat, jat) can be added to the cell_vector list.
                    for (const auto &it : distall[iat][jat]) {
                        if (exist[it.cell]) {
                            if (rc_tmp < 0.0 || it.dist <= rc_tmp) {
                                cell_vector.push_back(it.cell);
                            }
                        }
                    }
                    pairs_icell.push_back(cell_vector);
                }

                accum_tmp.clear();
                comb_cell.clear();
                cell_combination(pairs_icell, 0, accum_tmp, comb_cell);

                distance_list.clear();
                for (j = 0; j < comb_cell.size(); ++j) {

                    cellpair.clear();

                    for (k = 0; k < group_atom.size(); ++k) {
                        for (ii = 0; ii < group_atom[k]; ++ii) {
                            cellpair.push_back(comb_cell[j][k]);
                        }
                    }

                    dist_vector.clear();

                    for (k = 0; k < cellpair.size(); ++k) {
                        dist_tmp = distance(x_image[cellpair[k]][data_now[k]], x_image[0][iat]);
                        dist_vector.push_back(dist_tmp);
                    }

                    // Flag to check if the distance is smaller than the cutoff radius
                    isok = true;

                    for (k = 0; k < cellpair.size(); ++k) {
                        for (ii = k + 1; ii < cellpair.size(); ++ii) {
                            dist_tmp = distance(x_image[cellpair[k]][data_now[k]],
                                                x_image[cellpair[ii]][data_now[ii]]);
                            rc_tmp = cutoff_radii[order][kd[data_now[k]] - 1][kd[data_now[ii]] - 1];
                            if (rc_tmp >= 0.0 && dist_tmp > rc_tmp) {
                                isok = false;
                                break;
                            }
                            dist_vector.push_back(dist_tmp);
                        }
                        if (!isok) break;
                    }
                    if (isok) {
                        // This combination is a candidate of the minimum distance cluster
                        distance_list.emplace_back(MinDistList(cellpair, dist_vector));
                    }
                } // close loop over the mirror image combination

                if (!distance_list.empty()) {
                    // If the distance_list is not empty, there is a set of mirror images
                    // that satisfies the condition of the interaction.

                    pairs_icell.clear();
                    for (j = 0; j < intpair_uniq.size(); ++j) {
                        jat = intpair_uniq[j];
                        cell_vector.clear();

                        for (ii = 0; ii < mindist_pairs[iat][jat].size(); ++ii) {
                            cell_vector.push_back(mindist_pairs[iat][jat][ii].cell);
                        }
                        pairs_icell.push_back(cell_vector);
                    }

                    accum_tmp.clear();
                    comb_cell.clear();
                    comb_cell_atom_center.clear();
                    cell_combination(pairs_icell, 0, accum_tmp, comb_cell);

                    for (j = 0; j < comb_cell.size(); ++j) {
                        cellpair.clear();
                        for (k = 0; k < group_atom.size(); ++k) {
                            for (ii = 0; ii < group_atom[k]; ++ii) {
                                cellpair.push_back(comb_cell[j][k]);
                            }
                        }
                        comb_cell_atom_center.push_back(cellpair);
                    }

                    std::sort(distance_list.begin(), distance_list.end(),
                              MinDistList::compare_max_distance);
                    distmax = *std::max_element(distance_list[0].dist.begin(),
                                                distance_list[0].dist.end());
                    interaction_cluster_out[i].insert(InteractionCluster(data_now,
                                                                         comb_cell_atom_center,
                                                                         distmax));
                }
            }
        }
    }
    deallocate(list_now);
}


void Interaction::cell_combination(const std::vector<std::vector<int>> &array,
                                   const int i,
                                   std::vector<int> accum,
                                   std::vector<std::vector<int>> &comb)
{
    if (i == array.size()) {
        comb.push_back(accum);
    } else {
        std::vector<int> row = array[i];
        for (int j : row) {
            std::vector<int> tmp(accum);
            tmp.push_back(j);
            cell_combination(array, i + 1, tmp, comb);
        }
    }
}
