/*
 fcs.h

 Copyright (c) 2014--2017 Terumasa Tadano


 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm.h"
#include "interaction.h"
#include <vector>
#include <set>
#include <algorithm>

namespace ALM_NS
{
    class FcProperty
    {
    public:
        std::vector<int> elems; // flattened index of (iatom, icoordinate) in the supercell
        double sign; // factor (+1 or -1) to convert the mother FC to the child
        int mother; // index of the reducible force constants

        FcProperty();

        FcProperty(const FcProperty &obj) :
            elems(obj.elems), sign(obj.sign), mother(obj.mother) { }

        FcProperty(const int n,
                   const double c,
                   const int *arr,
                   const int m)
        {
            sign = c;
            mother = m;
            for (int i = 0; i < n; ++i) {
                elems.push_back(arr[i]);
            }
        }

        bool operator<(const FcProperty &a) const
        {
            return std::lexicographical_compare(elems.begin(), elems.end(),
                                                a.elems.begin(), a.elems.end());
        }

        bool operator==(const FcProperty &a) const
        {
            int n = elems.size();
            int n_ = a.elems.size();
            if (n != n_) return false;
            for (int i = 0; i < n; ++i) {
                if (elems[i] != a.elems[i]) return false;
            }
            return true;
        }
    };

    class ForceConstantTable
    {
    public:
        double fc_value;
        int multiplicity;
        std::vector<FcProperty> fclist;
        ForceConstantTable();
    };

    class Fcs
    {
    public:
        Fcs();
        ~Fcs();

        void init(ALM *);

        std::vector<int> *nequiv; // stores duplicate number of irreducible force constants
        std::vector<FcProperty> *fc_table; // all force constants
        std::vector<FcProperty> *fc_zeros; // zero force constants (due to space group symm.)

        void get_xyzcomponent(int, int **);

        bool is_inprim(int,
                       const int *,
                       int,
                       int **);

        int get_minimum_index_in_primitive(int,
                                           const int *,
                                           int,
                                           int,
                                           int **);
        double coef_sym(int,
                        double **,
                        const int *,
                        const int *);

        void generate_force_constant_table(int,
                                           int,
                                           const std::set<IntList> &,
                                           Symmetry *,
                                           std::string,
                                           std::vector<FcProperty> &,
                                           std::vector<int> &,
                                           std::vector<FcProperty> &,
                                           bool);

        void get_constraint_symmetry(int,
                                     Symmetry *,
                                     int,
                                     const std::set<IntList> &,
                                     std::string,
                                     const std::vector<FcProperty> &,
                                     int,
                                     double,
                                     std::vector<std::vector<double>> &);

    private:
        bool store_zeros;
        void set_default_variables();
        void deallocate_variables();
        bool is_ascending(int, const int *);
        bool is_inprim(int,
                       int,
                       int **);
        bool is_allzero(const std::vector<double> &,
                        double,
                        int &);
        void get_available_symmop(int,
                                  Symmetry *,
                                  std::string,
                                  int &,
                                  int **,
                                  double ***,
                                  bool);
    };
}

// Define a hash function for FcProperty class
// Use boost::hash_combine 
namespace std
{
    template <>
    struct hash<ALM_NS::FcProperty>
    {
        std::size_t operator ()(ALM_NS::FcProperty const &obj) const
        {
            hash<int> hasher;
            size_t seed = 0;
            for (auto i : obj.elems) {
                seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}
