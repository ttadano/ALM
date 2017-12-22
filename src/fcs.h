/*
 fcs.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>

namespace ALM_NS
{
    class FcProperty
    {
    public:
        std::vector<int> elems;
        double coef; // factor (usually +1 or -1)
        int mother; // index of the reducible force constants

        FcProperty();

        FcProperty(const FcProperty &obj) :
            coef(obj.coef), mother(obj.mother), elems(obj.elems)
        {
        }

        FcProperty(const int n, const double c, const int *arr, const int m)
        {
            coef = c;
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
    };

    class Fcs: protected Pointers
    {
    public:
        Fcs(class ALMCore *);
        ~Fcs();

        void init();

        std::vector<int> *ndup; // stores duplicate number of irreducible force constants
        std::vector<FcProperty> *fc_table; // all force constants

        std::string easyvizint(const int);
        void get_xyzcomponent(int, int **);
        void sort_tail(const int, int *);

        bool is_inprim(const int, const int *);
        bool is_inprim(const int);
        int min_inprim(const int, const int *);
        double coef_sym(const int, const int, const int *, const int *);

    private:
        void set_default_variables();
        void deallocate_variables();
        void generate_fclists(int maxorder, int *nzero);
        bool is_ascending(const int, const int *);
    };
}
