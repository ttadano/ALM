/*
 interaction.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
#include "constants.h"
#include "system.h"
#include "symmetry.h"
#include "timer.h"

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

        IntList(const IntList &a) : iarray(a.iarray) {};

        IntList(const int n,
                const int *arr)
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

        bool operator==(const IntList &a) const
        {
            int n = iarray.size();
            int n_ = a.iarray.size();
            if (n != n_) return false;
            for (int i = 0; i < n; ++i) {
                if (iarray[i] != a.iarray[i]) return false;
            }
            return true;
        }
    };


    class DistInfo
    {
    public:
        int cell;
        double dist;
        double relvec[3];

        DistInfo();

        DistInfo(const int n,
                 const double d,
                 const double x[3])
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

        DistList(const int n,
                 const double dist_tmp) : atom(n), dist(dist_tmp) { };

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
            : cell(cell_in), dist(dist_in) {};

        static bool compare_sum_distance(const MinDistList &a,
                                         const MinDistList &b)
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

        static bool compare_max_distance(const MinDistList &a,
                                         const MinDistList &b)
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

        InteractionCluster(const std::vector<int> &atom_in,
                           const std::vector<std::vector<int>> &cell_in,
                           const double dist_in)
            : atom(atom_in), cell(cell_in), distmax(dist_in) {};

        InteractionCluster(const std::vector<int> &atom_in,
                           const std::vector<std::vector<int>> &cell_in)
            : atom(atom_in), cell(cell_in), distmax(0.0) {};


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

        std::vector<std::string> str_order;
        std::set<IntList> *cluster_list;
        std::vector<int> **interaction_pair; // List of atoms inside the cutoff radius for each order
        std::set<InteractionCluster> **interaction_cluster;
        // Interaction many-body clusters with mirrow image information

        void init(const System *system,
                  const Symmetry *symmetry,
                  const int verbosity,
                  Timer *timer);

        bool satisfy_nbody_rule(const int,
                                const int *,
                                const int) const;

        bool is_incutoff(const int,
                         const int *,
                         const int,
                         const std::vector<int> &) const;

        void generate_interaction_information_by_cutoff(const int,
                                                        const int,
                                                        const std::vector<int> &,
                                                        const int * const *,
                                                        const double * const *,
                                                        std::vector<int> *) const;

        void set_interaction_by_cutoff(const unsigned int,
                                       const std::vector<int> &,
                                       const unsigned int,
                                       const int * const *,
                                       const double * const * const *,
                                       std::vector<int> **) const;
        void define(const int,
                    const unsigned int,
                    const int *,
                    const double * const * const *);
        int get_maxorder() const;
        int * get_nbody_include() const;

    private:

        int maxorder;
        int *nbody_include;
        double ***cutoff_radii;

        std::vector<DistInfo> **distall; // Distance of all pairs (i,j) under the PBC
        std::vector<DistInfo> **mindist_pairs; // All pairs (i,j) with the minimum distance

        void set_default_variables();
        void deallocate_variables();

        // can be made const function, but mindist_pairs is modified
        // in this function.
        void get_pairs_of_minimum_distance(const int,
                                           const double * const * const *,
                                           const int *);

        void print_neighborlist(const int,
                                const int,
                                const int * const *,
                                const std::vector<int> &,
                                const std::string *) const;

        void print_interaction_information(const int,
                                           const int * const *,
                                           const std::vector<int> &,
                                           const std::string *,
                                           const std::vector<int> * const *);

        void set_ordername();
        double distance(const double *, const double *) const;
        int nbody(const int, const int *) const;

        void calc_interaction_clusters(const int,
                                       const std::vector<int> &,
                                       const int * const *,
                                       const std::vector<int> * const *,
                                       const double * const *const *,
                                       const int *,
                                       std::set<InteractionCluster> **);

        void set_interaction_cluster(const int,
                                     const int,
                                     const std::vector<int> &,
                                     const int * const *,
                                     const std::vector<int> *,
                                     const double * const * const*,
                                     const int *,
                                     std::set<InteractionCluster> *) const;

        void cell_combination(const std::vector<std::vector<int>> &,
                              const int,
                              const std::vector<int>,
                              std::vector<std::vector<int>> &) const;

        void generate_pairs(const int,
                            const int * const *,
                            std::set<IntList> *,
                            const std::set<InteractionCluster> * const *) const;
    };
}

namespace std
{
    template <>
    struct hash<ALM_NS::IntList>
    {
        std::size_t operator ()(ALM_NS::IntList const &obj) const
        {
            hash<int> hasher;
            size_t seed = 0;
            for (auto i : obj.iarray) {
                seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}
