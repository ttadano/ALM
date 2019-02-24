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
        size_t number_of_atoms;
        size_t number_of_elems;
        std::vector<int> kind;
        std::vector<std::vector<double>> x_fractional;
        std::vector<std::vector<double>> x_cartesian;
    };

    class Spin
    {
    public:
        bool lspin = false;
        int time_reversal_symm = 1;
        int noncollinear = 0;
        std::vector<std::vector<double>> magmom;
    };

    class System
    {
    public:
        System();
        ~System();
        void init(const int,
                  Timer *);
        void frac2cart(double **) const;

        void set_supercell(const double [3][3],
                           const size_t,
                           const int *,
                           const double [][3],
                           const double transformation_matrix[3][3],
                           double primitive_axes[3][3]);
        void set_kdname(const std::string *);
        void set_periodicity(const int [3]);
        void set_spin_variables(const size_t nat,
                                const bool,
                                const int,
                                const int,
                                const double (*)[3]);
        void set_str_magmom(std::string);

        const Cell& get_supercell() const;
        double*** get_x_image() const;
        int* get_exist_image() const;
        std::string* get_kdname() const;
        int* get_periodicity() const;
        const Spin& get_spin() const;
        const std::string& get_str_magmom() const;
        const std::vector<std::vector<unsigned int>>& get_atomtype_group() const;

    private:
        // Variables for geometric structure
        Cell inputcell;
        Cell supercell, primcell;
        std::string *kdname;
        int *is_periodic; // is_periodic[3];
        double ***x_image;
        int *exist_image;

        // Variables for spins
        Spin inputspin, spin;
        std::string str_magmom;

        // concatenate atomic kind and magmom (only for collinear case)
        std::vector<std::vector<unsigned int>> atomtype_group;

        enum LatticeType { Direct, Reciprocal };

        void set_reciprocal_latt(const double [3][3],
                                 double [3][3]) const;
        void set_default_variables();
        void deallocate_variables();

        void set_inputcell(const double [3][3],
                           const size_t,
                           const int *,
                           const double [][3]);

        Cell generate_supercell(const Cell &cell_in,
                                const double transformation_matrix[3][3]);
                              
        void build_primitivecell(double primitive_axes[3][3],
                                 const bool refine_lattice_spglib,
                                 const double symprec_spglib);
        void find_primitive_spglib(const Cell &cell_in, Cell &primcell, 
                                   const bool refine_lattice,
                                   double primitive_axes[3][3],
                                   const double symprec);
        void get_transform_matrix_to_primitive(const std::string &symbol_in,
                                               double mat_out[3][3]);

        double volume(const double [3][3],
                      LatticeType) const;
        void set_atomtype_group();

        void generate_coordinate_of_periodic_images();
        void print_structure_stdout(const Cell &);
        void print_magmom_stdout() const;
    };
}
