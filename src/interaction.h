/*
 interaction.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm.h"
#include "constants.h"
#include <string>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>


namespace ALM_NS
{
    class IntList
    {
    public:
        std::vector<int> iarray;

        IntList()
        {
            iarray.clear();
        }

        IntList(const IntList &a)
        {
            for (auto p = a.iarray.begin(); p != a.iarray.end(); ++p) {
                iarray.push_back(*p);
            }
        }

        IntList(const int n, const int *arr)
        {
            for (int i = 0; i < n; ++i) {
                iarray.push_back(arr[i]);
            }
        }

        bool operator<(const IntList &a) const
        {
            return std::lexicographical_compare(iarray.begin(), iarray.end(),
                                                a.iarray.begin(), a.iarray.end());
        }
    };


    class DistInfo
    {
    public:
        int cell;
        double dist;
        double relvec[3];

        DistInfo();

        DistInfo(const int n, const double d, const double x[3])
        {
            cell = n;
            dist = d;
            for (int i = 0; i < 3; ++i) relvec[i] = x[i];
        }

        DistInfo(const DistInfo &obj)
        {
            cell = obj.cell;
            dist = obj.dist;
            for (int i = 0; i < 3; ++i) relvec[i] = obj.relvec[i];
        }

        bool operator<(const DistInfo &a) const
        {
            return dist < a.dist;
        }
    };

    class DistList
        // This class is used only in print_neighborlist. Can be replaced by a more generalic function.
    {
    public:
        int atom;
        double dist;

        DistList();

        DistList(const int n, const double dist_tmp) : atom(n), dist(dist_tmp)
        {
        };

        bool operator<(const DistList &a) const
        {
            if (std::abs(dist - a.dist) > eps8) {
                return dist < a.dist;
            }
            return atom < a.atom;
        }
    };

    class MinDistList
    {
    public:
        std::vector<int> cell;
        std::vector<double> dist;

        MinDistList();

        MinDistList(const std::vector<int> &cell_in,
                    const std::vector<double> &dist_in)
        {
            for (auto it = cell_in.begin(); it != cell_in.end(); ++it) {
                cell.push_back(*it);
            }
            for (auto it = dist_in.begin(); it != dist_in.end(); ++it) {
                dist.push_back(*it);
            }
        }

        static bool compare_sum_distance(const MinDistList &a, const MinDistList &b)
        {
            double dist_a = 0;
            double dist_b = 0;

            for (int i = 0; i < a.dist.size(); ++i) {
                dist_a += a.dist[i];
            }
            for (int i = 0; i < b.dist.size(); ++i) {
                dist_b += b.dist[i];
            }
            return dist_a < dist_b;
        }

        static bool compare_max_distance(const MinDistList &a, const MinDistList &b)
        {
            // This function works properly when dvec_a.size() > 0 and dvec_b.size() > 0
            std::vector<double> dvec_a, dvec_b;
            std::copy(a.dist.begin(), a.dist.end(), std::back_inserter(dvec_a));
            std::copy(b.dist.begin(), b.dist.end(), std::back_inserter(dvec_b));
            double max_dist_a = *std::max_element(dvec_a.begin(), dvec_a.end());
            double max_dist_b = *std::max_element(dvec_b.begin(), dvec_b.end());

            return max_dist_a < max_dist_b;
        }
    };

    class InteractionCluster
    {
    public:
        std::vector<int> atom;
        std::vector<std::vector<int>> cell;
        double distmax;

        InteractionCluster();

        InteractionCluster(const std::vector<int> atom_in,
                           const std::vector<std::vector<int>> cell_in,
                           const double dist_in)
        {
            for (int i = 0; i < atom_in.size(); ++i) {
                atom.push_back(atom_in[i]);
            }
            for (int i = 0; i < cell_in.size(); ++i) {
                cell.push_back(cell_in[i]);
            }
            distmax = dist_in;
        }

        InteractionCluster(const std::vector<int> atom_in,
                           const std::vector<std::vector<int>> cell_in)
        {
            for (int i = 0; i < atom_in.size(); ++i) {
                atom.push_back(atom_in[i]);
            }
            for (int i = 0; i < cell_in.size(); ++i) {
                cell.push_back(cell_in[i]);
            }
        }

        bool operator<(const InteractionCluster &a) const
        {
            return lexicographical_compare(atom.begin(), atom.end(),
                                           a.atom.begin(), a.atom.end());
        }
    };

    class Interaction
    {
    public:
        Interaction();
        ~Interaction();

        int maxorder;
        int *nbody_include;
        double ***rcs;

        std::vector<std::string> str_order;
        std::set<IntList> *cluster_list;
        std::vector<int> **interaction_pair; // List of atoms inside the cutoff radius for each order
        std::set<InteractionCluster> **interaction_cluster; // Interaction many-body clusters with mirrow image information

        void init(ALM *);

        bool satisfy_nbody_rule(const int,
                                const int *,
                                const int);

        bool is_incutoff(const int,
                         int *,
                         const int,
                         const std::vector<int>);

        void generate_interaction_information_by_cutoff(const int,
                                                        const int,
                                                        const std::vector<int> &,
                                                        int **,
                                                        double **,
                                                        std::vector<int> *);

        void set_interaction_by_cutoff(const unsigned int,
                                       const std::vector<int> &,
                                       const unsigned int,
                                       int **,
                                       double ***,
                                       std::vector<int> **);


    private:

        std::vector<DistInfo> **distall; // Distance of all pairs (i,j) under the PBC
        std::vector<DistInfo> **mindist_pairs; // All pairs (i,j) with the minimum distance

        void set_default_variables();
        void deallocate_variables();

        void get_pairs_of_minimum_distance(int,
                                           double ***,
                                           int *);

        void print_neighborlist(const int,
                                const int,
                                int **,
                                const std::vector<int> &,
                                std::string *);

        void print_interaction_information(const int,
                                           int **,
                                           const std::vector<int> &,
                                           std::string *,
                                           std::vector<int> **);

        void set_ordername();
        double distance(double *, double *);
        int nbody(const int, const int *);

        void calc_interaction_clusters(const int,
                                       const std::vector<int> &,
                                       int **,
                                       std::vector<int> **,
                                       double ***,
                                       int *,
                                       std::set<InteractionCluster> **);

        void set_interaction_cluster(const int,
                                     const int,
                                     const std::vector<int> &,
                                     int **,
                                     std::vector<int> *,
                                     double ***,
                                     int *,
                                     std::set<InteractionCluster> *);

        void cell_combination(std::vector<std::vector<int>>,
                              int, std::vector<int>,
                              std::vector<std::vector<int>> &);

        void generate_pairs(const int,
                            int **,
                            std::set<IntList> *,
                            std::set<InteractionCluster> **);
    };
}
