/*
cluster.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "cluster.h"
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

Cluster::Cluster()
{
    set_default_variables();
}

Cluster::~Cluster()
{
    deallocate_variables();
}

void Cluster::set_default_variables()
{
    maxorder = 0;
    nbody_include = nullptr;
    cutoff_radii = nullptr;
    distall = nullptr;
    mindist_pairs = nullptr;
    cluster_list = nullptr;
    interaction_pair = nullptr;
    cluster_each_atom = nullptr;
    cluster_table = nullptr;
}

void Cluster::deallocate_variables()
{
    if (nbody_include) {
        deallocate(nbody_include);
        nbody_include = nullptr;
    }
    if (cutoff_radii) {
        deallocate(cutoff_radii);
        cutoff_radii = nullptr;
    }
    if (cluster_list) {
        deallocate(cluster_list);
        cluster_list = nullptr;
    }
    if (mindist_pairs) {
        deallocate(mindist_pairs);
        mindist_pairs = nullptr;
    }
    if (interaction_pair) {
        deallocate(interaction_pair);
        interaction_pair = nullptr;
    }
    if (cluster_each_atom) {
        deallocate(cluster_each_atom);
        cluster_each_atom = nullptr;
    }
    if (distall) {
        deallocate(distall);
        distall = nullptr;
    }
    if (cluster_table) {
        deallocate(cluster_table);
        cluster_table = nullptr;
    }
}


void Cluster::init(const System *system,
                   const Symmetry *symmetry,
                   const int verbosity,
                   Timer *timer)
{
    timer->start_clock("cluster");

    int i;
    size_t j, k;
    const auto nat = system->get_supercell().number_of_atoms;
    const auto nkd = system->get_supercell().number_of_elems;

    if (verbosity > 0) {
        std::cout << " INTERACTION" << std::endl;
        std::cout << " ===========" << std::endl << std::endl;
    }

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
    allocate(interaction_pair, maxorder, symmetry->get_nat_prim());

    if (cluster_each_atom) {
        deallocate(cluster_each_atom);
    }
    allocate(cluster_each_atom, maxorder, symmetry->get_nat_prim());

    if (cluster_list) {
        deallocate(cluster_list);
    }
    allocate(cluster_list, maxorder);

    if (cluster_table) {
        deallocate(cluster_table);
    }
    allocate(cluster_table, maxorder);

    // If cutoff_radii and nbody_include is not set up to now,
    // use the default values.
    if (!cutoff_radii) {
        allocate(cutoff_radii, maxorder, nkd, nkd);
        for (i = 0; i < maxorder; ++i) {
            for (j = 0; j < nkd; ++j) {
                for (k = 0; k < nkd; ++k) {
                    cutoff_radii[i][j][k] = -1.0;
                }
            }
        }
    }
    if (!nbody_include) {
        allocate(nbody_include, maxorder);
        for (i = 0; i < maxorder; ++i) {
            nbody_include[i] = i + 2;
        }
    }

    generate_pairs_of_minimum_distance(nat,
                                       system->get_x_image(),
                                       system->get_exist_image());

    set_interaction_by_cutoff(system->get_supercell().number_of_atoms,
                              system->get_supercell().kind,
                              symmetry->get_nat_prim(),
                              symmetry->get_map_p2s());

    calc_interaction_clusters(symmetry->get_nat_prim(),
                              system->get_supercell().kind,
                              symmetry->get_map_p2s(),
                              system->get_x_image(),
                              system->get_exist_image());

    identify_irreducible_clusters(symmetry, symmetry->get_map_sym());

    if (verbosity > 0) {
        std::cout << "  +++ Cutoff Radii Matrix (NKD x NKD matrix) +++" << std::endl;

        for (i = 0; i < maxorder; ++i) {
            std::cout << "  " << std::setw(9) << get_ordername(i) << std::endl;
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

        print_neighborlist(system->get_supercell().number_of_atoms,
                           symmetry->get_nat_prim(),
                           symmetry->get_map_p2s(),
                           system->get_supercell().kind,
                           system->get_kdname());
    }

    if (verbosity > 1) {
        print_cluster_information(symmetry->get_nat_prim(),
                                  symmetry->get_map_p2s(),
                                  system->get_supercell().kind,
                                  system->get_kdname(),
                                  system->get_x_image());
    }

    if (verbosity > 0) {
        for (i = 0; i < maxorder; ++i) {
            if (i + 2 > nbody_include[i]) {
                std::cout << "  For " << std::setw(8) << get_ordername(i) << ", ";
                std::cout << "interactions related to more than" << std::setw(2) << nbody_include[i];
                std::cout << " atoms will be neglected." << std::endl;
            }
        }


        timer->print_elapsed();
        std::cout << " -------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
    }

    timer->stop_clock("cluster");
}


void Cluster::generate_pairs_of_minimum_distance(const size_t nat,
                                                 const double *const *const *xc_in,
                                                 const int *exist) const
{
    size_t i, j;
    double vec[3];

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < nat; ++j) {

            distall[i][j].clear();

            for (auto icell = 0; icell < 27; ++icell) {

                if (exist[icell]) {

                    const auto dist_tmp = distance(xc_in[0][i], xc_in[icell][j]);

                    for (auto k = 0; k < 3; ++k) vec[k] = xc_in[icell][j][k] - xc_in[0][i][k];

                    distall[i][j].emplace_back(DistInfo(icell, dist_tmp, vec));
                }
            }
            std::sort(distall[i][j].begin(), distall[i][j].end());
        }
    }

    // Construct pairs of minimum distance.

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < nat; ++j) {
            mindist_pairs[i][j].clear();

            const auto dist_min = distall[i][j][0].dist;
            for (const auto &it : distall[i][j]) {
                // The tolerance below (1.e-3) should be chosen so that
                // the mirror images with equal distances are found correctly.
                // If this fails, the phonon dispersion would be incorrect.
                if (std::abs(it.dist - dist_min) < 1.0e-3) {
                    mindist_pairs[i][j].emplace_back(it);
                }
            }
        }
    }
}

void Cluster::set_interaction_by_cutoff(const size_t nat,
                                        const std::vector<int> &kd,
                                        const size_t nat_prim,
                                        const std::vector<std::vector<int>> &map_p2s) const
{
    // For each atom 'i' in the unit cell at origin, generate a list of atoms 'j'
    // whose distance from atom 'i' is smaller than the given cutoff radius.

    for (auto order = 0; order < maxorder; ++order) {

        for (size_t i = 0; i < nat_prim; ++i) {

            interaction_pair[order][i].clear();

            const auto iat = map_p2s[i][0];
            const auto ikd = kd[iat] - 1;

            for (size_t jat = 0; jat < nat; ++jat) {

                const auto jkd = kd[jat] - 1;
                const auto cutoff_tmp = cutoff_radii[order][ikd][jkd];

                if (cutoff_tmp < 0.0) {

                    // Cutoff 'None'
                    interaction_pair[order][i].push_back(jat);

                } else {

                    if (mindist_pairs[iat][jat][0].dist <= cutoff_tmp) {
                        interaction_pair[order][i].push_back(jat);
                    }
                }
            }
        }
    }
}

void Cluster::print_neighborlist(const size_t nat,
                                 const size_t natmin,
                                 const std::vector<std::vector<int>> &map_p2s,
                                 const std::vector<int> &kd,
                                 const std::string *kdname) const
{
    //
    // Print the list of neighboring atoms and distances
    //
    size_t j, k;
    int icount;

    std::vector<std::pair<size_t, double>> neighbor;
    std::pair<size_t, double> pair_tmp;
    std::vector<int> atomlist;

    std::cout << std::endl;
    std::cout << "  List of neighboring atoms below." << std::endl;
    std::cout << "  Format [N th-nearest shell, distance (Number of atoms on the shell)]"
        << std::endl << std::endl;

    for (size_t i = 0; i < natmin; ++i) {

        auto iat = map_p2s[i][0];
        neighbor.clear();

        for (j = 0; j < nat; ++j) {
            pair_tmp.first = j;
            pair_tmp.second = mindist_pairs[iat][j][0].dist;
            neighbor.push_back(pair_tmp);
        }
        std::sort(neighbor.begin(), neighbor.end(),
                  [](const std::pair<size_t, double> &a,
                     const std::pair<size_t, double> &b)
                  {
                      if (std::abs(a.second - b.second) > eps8) {
                          return a.second < b.second;
                      }
                      return a.first < b.first;
                  });


        auto nthnearest = 0;
        atomlist.clear();

        iat = map_p2s[i][0];
        std::cout << std::setw(5) << iat + 1 << " ("
            << std::setw(3) << kdname[kd[iat] - 1] << "): ";

        auto dist_tmp = 0.0;

        for (j = 0; j < nat; ++j) {

            if (neighbor[j].second < eps8) continue; // distance is zero

            if (std::abs(neighbor[j].second - dist_tmp) > eps6) {

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


                dist_tmp = neighbor[j].second;
                atomlist.clear();
                atomlist.push_back(neighbor[j].first);
            } else {
                atomlist.push_back(neighbor[j].first);
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
}


int Cluster::get_maxorder() const
{
    return maxorder;
}

void Cluster::define(const int maxorder_in,
                     const size_t nkd,
                     const int *nbody_include_in,
                     const double *cutoff_radii_in)
{
    maxorder = maxorder_in;
    if (nbody_include) {
        deallocate(nbody_include);
    }
    allocate(nbody_include, maxorder);

    for (auto i = 0; i < maxorder; ++i) {
        nbody_include[i] = nbody_include_in[i];
    }

    // nkd == 0: Use current values, which are probably default values
    // defined at Cluster::init().
    if (nkd > 0) {
        if (cutoff_radii) {
            deallocate(cutoff_radii);
        }
        allocate(cutoff_radii, maxorder, nkd, nkd);
    }

    // if cutoff_radii_in = nullptr, use current value
    if (cutoff_radii_in) {
        auto counter = 0;
        for (auto i = 0; i < maxorder; ++i) {
            for (size_t j = 0; j < nkd; ++j) {
                for (size_t k = 0; k < nkd; ++k) {
                    cutoff_radii[i][j][k] = cutoff_radii_in[counter++];
                }
            }
        }
    }
}

int* Cluster::get_nbody_include() const
{
    return nbody_include;
}

std::string Cluster::get_ordername(const unsigned int order) const
{
    if (order == 0) {
        return "HARMONIC";
    } else {
        return "ANHARM" + std::to_string(order + 2);
    }
}

const std::set<ClusterEntry>& Cluster::get_cluster_list(const unsigned int order) const
{
    return cluster_list[order];
}

const std::vector<int>& Cluster::get_interaction_pair(const unsigned int order,
                                                      const size_t atom_index) const
{
    return interaction_pair[order][atom_index];
}

const std::set<InteractionCluster>& Cluster::get_cluster_each_atom(const unsigned int order,
                                                                   const size_t atom_index) const
{
    return cluster_each_atom[order][atom_index];
}

std::vector<double> Cluster::get_cluster_distance(const std::vector<int> &atoms,
                                                  const std::vector<int> &kd,
                                                  const double *const*const*x_image,
                                                  const bool use_cutoff_radii) const
{
    std::vector<double> distances(2);
    auto atoms_unique(atoms);

    const auto order = atoms.size() - 2;

    // Remove duplicate entries
    std::sort(atoms_unique.begin(), atoms_unique.end());
    atoms_unique.erase(std::unique(atoms_unique.begin(), atoms_unique.end()),
                       atoms_unique.end());

    if (atoms_unique.size() == 1) {
        distances[0] = 0.0;
        distances[1] = 0.0;
        return distances;
    }

    // For multi-body terms, searching the maximum and minimum distances of a cluster
    // under the periodic boundary condition is not trivial.

    std::vector<int> cell_vector, accum_tmp;
    std::vector<std::vector<int>> pairs_icell, comb_cell;

    const auto iat = atoms_unique[0];

    for (size_t i = 0; i < atoms_unique.size() - 1; ++i) {
        const auto jat = atoms_unique[i + 1];

        cell_vector.clear();

        // Loop over the cell images of atom 'jat' and add to the list
        // as a candidate for the cluster.

        if (use_cutoff_radii) {
            // The mirror images whose distance is larger than the minimum value
            // of the distance(iat, jat) may also be added to the cell_vector list.
            const auto rc_tmp = cutoff_radii[order][kd[iat] - 1][kd[jat] - 1];
            for (const auto &it : distall[iat][jat]) {
                if (rc_tmp < 0.0 || it.dist <= rc_tmp) {
                    cell_vector.push_back(it.cell);
                }
            }
        } else {
            for (const auto &it : distall[iat][jat]) {
                cell_vector.push_back(it.cell);
            }
        }

        pairs_icell.push_back(cell_vector);
    }

    accum_tmp.clear();
    comb_cell.clear();
    cell_combination(pairs_icell, 0, accum_tmp, comb_cell);

    std::vector<double> dist_vector;

    auto dist_max = 1.0e+8;
    auto dist_min = 1.0e+8;

    for (const auto &cell_now : comb_cell) {

        dist_vector.clear();

        for (size_t j = 0; j < cell_now.size(); ++j) {
            const auto dist_tmp = distance(x_image[cell_now[j]][atoms_unique[j + 1]],
                                           x_image[0][iat]);
            dist_vector.push_back(dist_tmp);
        }

        for (size_t j = 0; j < cell_now.size(); ++j) {
            for (size_t k = j + 1; k < cell_now.size(); ++k) {
                const auto dist_tmp = distance(x_image[cell_now[j]][atoms_unique[j + 1]],
                                               x_image[cell_now[k]][atoms_unique[k + 1]]);
                dist_vector.push_back(dist_tmp);
            }
        }
        std::sort(dist_vector.begin(), dist_vector.end());

        // If the maximum distance is smaller than the previously found one,
        // the new one can be a candidate of the combination of mirror images
        // that gives the minimum cluster size.
        if (dist_vector.back() < dist_max) {
            dist_max = dist_vector.back();

            // If the minimum distance is smaller as well,
            // it should be the smaller cluster size.
            if (dist_vector.front() < dist_min) {
                dist_min = dist_vector.front();
            }
        }
    }
    distances[0] = dist_min;
    distances[1] = dist_max;

    return distances;
}

//void Cluster::erase_cluster(const unsigned order,
//                            const std::vector<int> &target_cluster,
//                            const Symmetry *symmetry)
//{
//    if (order >= maxorder) return;
//
//    for (auto i = 0; i < target_cluster.size(); ++i) {
//
//        auto iat = target_cluster[i];
//        auto is_inside_primitive = false;
//
//        for (auto j = 0; j < symmetry->get_nat_prim(), ++j) {
//            
//        }
//    
//    }
//}


void Cluster::print_cluster_information(const size_t natmin,
                                        const std::vector<std::vector<int>> &map_p2s,
                                        const std::vector<int> &kd,
                                        const std::string *kdname,
                                        const double *const*const*x_image) const
{
    std::vector<int> intlist;

    std::cout << std::endl;
    std::cout << "  List of irreducible clusters:" << std::endl;


    for (auto order = 0; order < maxorder; ++order) {

        std::cout << std::endl << "   ***" << get_ordername(order) << "***\n\n";

        std::cout << "   Number of unique clusters:" << std::setw(6) << cluster_table[order].size() << std::endl;
        std::cout << "   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";


        std::vector<SortCluster> cluster_sorted;

        for (const auto &it : cluster_table[order]) {
            const auto it2 = it.begin();
            const auto distances = get_cluster_distance(it2->iarray, kd, x_image, true);
            cluster_sorted.push_back(SortCluster(it2->iarray,
                                                 it.size(),
                                                 distances[0],
                                                 distances[1]));
        }

        // Sort twice
        std::stable_sort(cluster_sorted.begin(), cluster_sorted.end(),
                         [](const SortCluster &a,
                            const SortCluster &b)
                         {
                             return a.distance_min < b.distance_min;
                         });
        std::stable_sort(cluster_sorted.begin(), cluster_sorted.end(),
                         [](const SortCluster &a,
                            const SortCluster &b)
                         {
                             return a.distance_max < b.distance_max;
                         });

        size_t i = 0;

        for (const auto &it : cluster_sorted) {


            std::cout << "   " << std::setw(6) << i + 1;
            std::cout << "(" << std::setw(5) << it.multiplicity << ") | ";

            for (const auto j : it.atoms) {
                std::cout << std::setw(5) << j + 1;
                std::cout << "(" << std::setw(3) << kdname[kd[j] - 1] << ")";
            }
            std::cout << std::setw(10) << it.distance_min;
            std::cout << std::setw(10) << it.distance_max;
            std::cout << "\n";

            ++i;
        }
        std::cout << "   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    }

    std::cout << std::endl;
}


bool Cluster::is_incutoff(const int n,
                          const int *atomnumlist,
                          const size_t order,
                          const std::vector<int> &kd) const
{
    for (auto i = 0; i < n; ++i) {
        const auto iat = atomnumlist[i];
        const auto ikd = kd[iat] - 1;

        for (auto j = i + 1; j < n; ++j) {
            const auto jat = atomnumlist[j];
            const auto jkd = kd[jat] - 1;

            const auto cutoff_tmp = cutoff_radii[order][ikd][jkd];

            if (cutoff_tmp >= 0.0 &&
                (mindist_pairs[iat][jat][0].dist > cutoff_tmp)) {
                return false;
            }

        }
    }
    return true;
}

int Cluster::nbody(const int n,
                   const int *arr) const
{
    std::vector<int> v(n);

    for (auto i = 0; i < n; ++i) {
        v[i] = arr[i];
    }
    std::stable_sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());

    const auto ret = v.size();
    v.clear();

    return static_cast<int>(ret);
}

bool Cluster::satisfy_nbody_rule(const int nelem,
                                 const int *arr,
                                 const int order) const
{
    return nbody(nelem, arr) <= nbody_include[order];
}


void Cluster::calc_interaction_clusters(const size_t natmin,
                                        const std::vector<int> &kd,
                                        const std::vector<std::vector<int>> &map_p2s,
                                        const double *const *const *x_image,
                                        const int *exist) const
{
    //
    // Calculate the complete set of clusters for all orders.
    //

    std::vector<std::vector<int>> **cluster_tmp;

    allocate(cluster_tmp, maxorder, natmin);
    for (auto order = 0; order < maxorder; ++order) {

        search_clusters(order,
                        natmin,
                        kd,
                        map_p2s,
                        interaction_pair[order],
                        x_image,
                        exist,
                        cluster_tmp[order]);

        // This method is left here temporary and will be removed in future.
        set_interaction_cluster(order,
                                natmin,
                                kd,
                                map_p2s,
                                interaction_pair[order],
                                x_image,
                                exist,
                                cluster_each_atom[order]);

    }
    deallocate(cluster_tmp);
}


void Cluster::set_interaction_cluster(const int order,
                                      const size_t natmin,
                                      const std::vector<int> &kd,
                                      const std::vector<std::vector<int>> &map_p2s,
                                      const std::vector<int> *interaction_pair_in,
                                      const double *const *const *x_image,
                                      const int *exist,
                                      std::set<InteractionCluster> *interaction_cluster_out) const
{
    //
    // Calculate a set of clusters for the given order
    //

    size_t j, k;
    int jat;
    size_t ii;
    int jkd;
    int icount;

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

    for (size_t i = 0; i < natmin; ++i) {

        interaction_cluster_out[i].clear();

        const auto iat = map_p2s[i][0];
        const auto ikd = kd[iat] - 1;
        list_now[0] = iat;

        // List of 2-body interaction pairs
        auto intlist(interaction_pair_in[i]);
        std::sort(intlist.begin(), intlist.end()); // Need to sort here

        if (order == 0) {

            // Harmonic term

            for (auto ielem : intlist) {

                jat = ielem;
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

            // First, we generate all possible combinations of clusters.
            CombinationWithRepetition<int> g(intlist.begin(), intlist.end(), order + 1);
            do {
                auto data = g.now();

                list_now[0] = iat;
                for (auto m = 0; m < order + 1; ++m) list_now[m + 1] = data[m];

                // Save as a candidate if the cluster satisfies the NBODY-rule.
                if (satisfy_nbody_rule(order + 2, list_now, order)) {
                    data_vec.emplace_back(data);
                }

            } while (g.next());

            intlist.clear();

            const auto ndata = data_vec.size();

            for (size_t idata = 0; idata < ndata; ++idata) {

                data_now = data_vec[idata];

                // Uniq the list of atoms in data like as follows:
                // cubic   term : (i, i) --> (i) x 2
                // quartic term : (i, i, j) --> (i, j) x (2, 1)
                intpair_uniq.clear();
                group_atom.clear();
                icount = 1;

                for (auto m = 0; m < order; ++m) {
                    if (data_now[m] == data_now[m + 1]) {
                        ++icount;
                    } else {
                        group_atom.push_back(icount);
                        intpair_uniq.push_back(data_now[m]);
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
                    // as a candidate for the cluster.
                    // The mirror images whose distance is larger than the minimum value
                    // of the distance(iat, jat) can be added to the cell_vector list.
                    for (const auto &it : distall[iat][jat]) {
                        if (rc_tmp < 0.0 || it.dist <= rc_tmp) {
                            cell_vector.push_back(it.cell);
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
                        for (auto m = 0; m < group_atom[k]; ++m) {
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
                    // that satisfies the condition of the cluster.

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
                            for (auto m = 0; m < group_atom[k]; ++m) {
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


void Cluster::cell_combination(const std::vector<std::vector<int>> &array,
                               const size_t i,
                               const std::vector<int> &accum,
                               std::vector<std::vector<int>> &comb)
{
    if (i == array.size()) {
        comb.push_back(accum);
    } else {
        auto row = array[i];
        for (auto j : row) {
            auto tmp(accum);
            tmp.push_back(j);
            cell_combination(array, i + 1, tmp, comb);
        }
    }
}


double Cluster::distance(const double *x1,
                         const double *x2) const
{
    auto dist = std::pow(x1[0] - x2[0], 2)
        + std::pow(x1[1] - x2[1], 2)
        + std::pow(x1[2] - x2[2], 2);
    dist = std::sqrt(dist);

    return dist;
}


void Cluster::search_clusters(const int order,
                              const size_t natmin,
                              const std::vector<int> &kd,
                              const std::vector<std::vector<int>> &map_p2s,
                              const std::vector<int> *interaction_pair_in,
                              const double *const*const*x_image,
                              const int *exist,
                              std::vector<std::vector<int>> *cluster_tmp) const
{
    //
    // Calculate a set of clusters for the given order
    //

    size_t j, k;
    int jat;
    size_t ii;
    int jkd;
    int icount;

    double dist_tmp, rc_tmp;

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
    // std::vector<MinDistList> distance_list;
    std::vector<std::vector<int>> data_vec;

    std::vector<ClusterEntry> cluster_merged;

    cluster_merged.clear();

    allocate(list_now, order + 2);

    for (size_t i = 0; i < natmin; ++i) {

        cluster_tmp[i].clear();

        const auto iat = map_p2s[i][0];
        const auto ikd = kd[iat] - 1;
        list_now[0] = iat;

        // List of 2-body interaction pairs
        auto intlist(interaction_pair_in[i]);
        std::sort(intlist.begin(), intlist.end()); // Need to sort here

        if (order == 0) {

            // Harmonic term

            for (auto ielem : intlist) {

                jat = ielem;
                list_now[1] = jat;

                if (!satisfy_nbody_rule(2, list_now, 0)) continue;

                comb_cell_min.clear();
                atom_tmp.clear();
                atom_tmp.push_back(jat);

                cluster_tmp[i].emplace_back(atom_tmp);
            }

        } else if (order > 0) {

            // Anharmonic terms

            data_vec.clear();

            // First, we generate all possible combinations of clusters.
            CombinationWithRepetition<int> g(intlist.begin(), intlist.end(), order + 1);
            do {
                auto data = g.now();

                list_now[0] = iat;
                for (auto m = 0; m < order + 1; ++m) list_now[m + 1] = data[m];

                // Save as a candidate if the cluster satisfies the NBODY-rule.
                if (satisfy_nbody_rule(order + 2, list_now, order)) {
                    data_vec.emplace_back(data);
                }

            } while (g.next());

            intlist.clear();

            const auto ndata = data_vec.size();

            for (size_t idata = 0; idata < ndata; ++idata) {

                data_now = data_vec[idata];

                // Uniq the list of atoms in data like as follows:
                // cubic   term : (i, i) --> [i] x 2
                // quartic term : (i, i, j) --> [i] x 2 + [j] x 1
                intpair_uniq.clear();
                group_atom.clear();
                icount = 1;

                for (auto m = 0; m < order; ++m) {
                    if (data_now[m] == data_now[m + 1]) {
                        ++icount;
                    } else {
                        group_atom.push_back(icount);
                        intpair_uniq.push_back(data_now[m]);
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
                    // as a candidate for the cluster.
                    // The mirror images whose distance is larger than the minimum value
                    // of the distance(iat, jat) may also be added to the cell_vector list.
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

                for (j = 0; j < comb_cell.size(); ++j) {

                    cellpair.clear();

                    for (k = 0; k < group_atom.size(); ++k) {
                        for (auto m = 0; m < group_atom[k]; ++m) {
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
                        // If any of the mirror image combination satisfies the 
                        // cutoff radii criterion, add the cluster to the list and break the loop.
                        cluster_tmp[i].emplace_back(data_now);
                        break;
                    }
                } // close loop over the mirror image combination
            }
        }

        for (const auto &it : cluster_tmp[i]) {

            list_now[0] = iat;
            for (j = 0; j < order + 1; ++j) {
                list_now[j + 1] = it[j];
            }
            insort(order + 2, list_now);
            cluster_merged.emplace_back(order + 2, list_now);
        }
    }

    std::sort(cluster_merged.begin(), cluster_merged.end());
    cluster_merged.erase(std::unique(cluster_merged.begin(), cluster_merged.end()),
                         cluster_merged.end());

    // Create cluster_list for later reference.
    cluster_list[order].clear();

    for (const auto &it : cluster_merged) {
        cluster_list[order].insert(it);
    }
    deallocate(list_now);
}

void Cluster::identify_irreducible_clusters(const Symmetry *symm_in,
                                            const std::vector<std::vector<int>> &map_sym) const
{
    const auto nsym = map_sym[0].size();

    std::vector<int> atoms, atoms_symm;
    size_t i;
    std::vector<ClusterEntry> cluster_table_tmp;

    for (size_t order = 0; order < maxorder; ++order) {
        atoms.resize(order + 2);
        atoms_symm.resize(order + 2);

        auto cluster_list_copy(cluster_list[order]);
        auto it_orig = cluster_list_copy.begin();

        cluster_table[order].clear();

        while (it_orig != cluster_list_copy.end()) {

            for (i = 0; i < order + 2; ++i) {
                atoms[i] = it_orig->iarray[i];
            }

            cluster_table_tmp.clear();

            for (size_t isym = 0; isym < nsym; ++isym) {
                for (i = 0; i < order + 2; ++i) {
                    atoms_symm[i] = map_sym[atoms[i]][isym];
                }

                if (!is_ascending(atoms_symm)) continue;
                if (!symm_in->is_anyof_inside_primitive(atoms_symm)) continue;

                auto it_target = cluster_list_copy.find(ClusterEntry(order + 2, &atoms_symm[0]));

                // If the newly generated cluster is a member of the original cluster,
                // remove it from the cluster_list_copy.
                if (it_target != cluster_list_copy.end()) {
                    it_orig = cluster_list_copy.erase(it_target);
                }

                // Add the current cluster to a temporary list
                cluster_table_tmp.push_back(ClusterEntry(order + 2, &atoms_symm[0]));
            }

            // Sort & uniq the temporary cluster
            std::sort(cluster_table_tmp.begin(), cluster_table_tmp.end());
            cluster_table_tmp.erase(std::unique(cluster_table_tmp.begin(), cluster_table_tmp.end()),
                                    cluster_table_tmp.end());

            // Register the temporary cluster to the cluster_table
            cluster_table[order].push_back(cluster_table_tmp);

        }
    }
}

bool Cluster::is_ascending(const std::vector<int> &arr)
{
    const auto n = arr.size();
    for (size_t i = 0; i < n - 1; ++i) {
        if (arr[i] > arr[i + 1]) return false;
    }
    return true;
}
