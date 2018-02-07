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

    class Cell
    {
    public:
        double lattice_vector[3][3];
        double reciprocal_lattice_vector[3][3];
        double volume;
        unsigned int number_of_atmos;
        unsigned int number_of_elems;
        std::vector<int> kind;
        std::vector<std::vector<double>> x_fractional;
        std::vector<std::vector<double>> x_cartesian;
    };


    class System
    {
    public:
        System();
        ~System();
        void init(ALM *);
        void frac2cart(double **);
        void load_reference_system_xml(Symmetry *, Fcs *, std::string, const int, double *);

        void set_cell(const double [3][3], const unsigned int,
                      const unsigned int, int *, double **, Cell &);

        Cell primitivecell, supercell;
        int nat, nkd;
        int ndata, nstart, nend;
        int *kd;
        double lavec[3][3];
        double **xcoord; // fractional coordinate

        std::string *kdname;
        unsigned int nclassatom;
        std::vector<unsigned int> *atomlist_class;

        bool lspin;
        int trev_sym_mag;
        double **magmom;
        std::string str_magmom;
        int noncollinear;

    private:
        enum LatticeType { Direct, Reciprocal };

        void set_reciprocal_latt(const double [3][3], double [3][3]);
        void set_default_variables();
        void deallocate_variables();

        double volume(const double [3][3], LatticeType);
        void setup_atomic_class(int *);
        void print_structure_stdout(const Cell &);
        void print_magmom_stdout();
    };
}
