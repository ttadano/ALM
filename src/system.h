/*
 system.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <vector>
#include "timer.h"

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

    class Spin
    {
    public:
        bool lspin;
        int time_reversal_symm;
        int noncollinear;
        std::vector<std::vector<double>> magmom;
    };

    class System
    {
    public:
        System();
        ~System();
        void init(const int verbosity,
                  Timer *timer);
        void frac2cart(double **) const;

        void set_supercell(const double [3][3],
                           const unsigned int,
                           const unsigned int,
                           const int *,
                           const double [][3]);

        Cell get_supercell() const;

        void set_spin_variable(const bool,
                               const int,
                               const int,
                               const unsigned int,
                               double **);

        Spin spin;

        std::string *kdname;
        // concatenate atomic kind and magmom (only for collinear case)
        std::vector<std::vector<unsigned int>> atomtype_group;
        int is_periodic[3];
        double ***x_image;
        int *exist_image;

        // Variables for spins

        bool lspin;
        int trev_sym_mag;
        int noncollinear;
        double **magmom;
        std::string str_magmom;

    private:
        Cell supercell;

        enum LatticeType { Direct, Reciprocal };

        void set_reciprocal_latt(const double [3][3], double [3][3]) const;
        void set_default_variables();
        void deallocate_variables();

        double volume(const double [3][3], LatticeType) const;
        void set_atomtype_group();

        void generate_coordinate_of_periodic_images(const unsigned int,
                                                    const std::vector<std::vector<double>> &,
                                                    const int [3],
                                                    double ***,
                                                    int *) const;
        void print_structure_stdout(const Cell &);
        void print_magmom_stdout() const;
    };
}
