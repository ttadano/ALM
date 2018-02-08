/*
 system.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm.h"
#include <string>
#include <vector>

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
            }
            if (this->element == a.element) {
                return this->magmom < a.magmom;
            }
            return false;
        }
    };

    class Cell
    {
    public:
        double lattice_vector[3][3];
        double reciprocal_lattice_vector[3][3];
        double volume;
        unsigned int number_of_atoms;
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

        void set_cell(const double [3][3], unsigned int,
                      unsigned int, int *, double **, Cell &);

        Cell primitivecell, supercell;

        std::string *kdname;
        unsigned int nclassatom;
        std::vector<unsigned int> *atomlist_class;
        int is_periodic[3];

        // Variables for spins
        bool lspin;
        int trev_sym_mag;
        int noncollinear;
        double **magmom;
        std::string str_magmom;

        // Referenced from input_setter, writer
        int nat, nkd;
        int *kd;
        double lavec[3][3];
        double **xcoord; // fractional coordinate

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
