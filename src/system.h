/*
 system.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

//#include "pointers.h"
#include <string>
#include <vector>
#include "alm.h"

namespace ALM_NS
{
    class AtomType
    {
    public:
        int element;
        double magmom;

        bool operator<(const AtomType &a) const
        {
            if (this->element < a.element) {
                return true;
            } else if (this->element == a.element) {
                return this->magmom < a.magmom;
            } else {
                return false;
            }
        }
    };


    class System
    {
    public:
        System();
        ~System();
        void init(ALM *);
        void recips(double [3][3], double [3][3]);
        void frac2cart(double **);
        void load_reference_system();
        void load_reference_system_xml(std::string, const int, double *);

        int nat, nkd;
        //        int nat_prim;
        int ndata, nstart, nend;
        int *kd;
        //        int *kd_prim;
        double lavec[3][3], rlavec[3][3];
        double **xcoord; // fractional coordinate
        double **x_cartesian;
        //        double lavec_prim[3][3], rlavec_prim[3][3];
        double **magmom;
        std::string str_magmom;
        int noncollinear;
        std::string *kdname;

        unsigned int nclassatom;

        std::vector<unsigned int> *atomlist_class;
        bool lspin;
        double cell_volume;


    private:
        void set_default_variables();
        void deallocate_variables();
        double volume(double [3], double [3], double [3]);
        void setup_atomic_class(int *);
    };
}
