/*
 constraint.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm.h"
#include "constants.h"
#include "fcs.h"
#include "system.h"
#include <boost/bimap.hpp>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

namespace ALM_NS
{
    class ConstraintClass
    {
    public:
        std::vector<double> w_const;

        ConstraintClass();

        ConstraintClass(const ConstraintClass &a) : w_const(a.w_const) { }

        ConstraintClass(std::vector<double> vec) : w_const(std::move(vec)) { }

        ConstraintClass(const int n,
                        const double *arr,
                        const int nshift = 0)
        {
            for (int i = nshift; i < n; ++i) {
                w_const.push_back(arr[i]);
            }
        }

        bool operator<(const ConstraintClass &a) const
        {
            return std::lexicographical_compare(w_const.begin(), w_const.end(),
                                                a.w_const.begin(), a.w_const.end());
        }
    };

    class ConstraintTypeFix
    {
    public:
        unsigned int p_index_target;
        double val_to_fix;

        ConstraintTypeFix(const unsigned int index_in,
                          const double val_in) :
            p_index_target(index_in), val_to_fix(val_in) { }
    };

    class ConstraintTypeRelate
    {
    public:
        unsigned int p_index_target;
        std::vector<double> alpha;
        std::vector<unsigned int> p_index_orig;

        ConstraintTypeRelate(const unsigned int index_in,
                             const std::vector<double> &alpha_in,
                             const std::vector<unsigned int> &p_index_in) :
            p_index_target(index_in), alpha(alpha_in), p_index_orig(p_index_in) { }
    };

    inline bool equal_within_eps12(const std::vector<double> &a,
                                   const std::vector<double> &b)
    {
        int n = a.size();
        int m = b.size();
        if (n != m) return false;
        double res = 0.0;
        for (int i = 0; i < n; ++i) {
            if (std::abs(a[i] - b[i]) > eps12) return false;
        }
        //        if (std::sqrt(res)>eps12) return false;
        return true;
    }

    class ConstraintIntegerElement
    {
        // For sparse representation
    public:
        unsigned int col;
        int val;

        ConstraintIntegerElement(const unsigned int col_in,
                                 const int val_in) : 
                                 col(col_in), val(val_in) {}
    };

    inline bool operator<(const std::vector<ConstraintIntegerElement> &obj1, 
                          const std::vector<ConstraintIntegerElement> &obj2) {
                
                const int len1 = obj1.size();
                const int len2 = obj2.size();
                const int min = std::min(len1, len2);

                for (int i = 0; i < min; ++i) {
                    if (obj1[i].col < obj2[i].col) {
                        return true;
                    } else if (obj1[i].col > obj2[i].col) {
                        return false;
                    } else {
                        if (obj1[i].val < obj2[i].val) {
                            return true;
                        } else if (obj1[i].val > obj2[i].val) {
                            return false;
                        }
                    }
                }
                return false;
    }

    inline bool operator==(const std::vector<ConstraintIntegerElement> &obj1, 
                        const std::vector<ConstraintIntegerElement> &obj2) {
            
            const int len1 = obj1.size();
            const int len2 = obj2.size();
            if (len1 != len2) return false;

            for (int i = 0; i < len1; ++i) {
                if (obj1[i].col != obj2[i].col || obj1[i].val != obj2[i].val) {
                    return false;
                } 
            }
            return true;
    }

    class Constraint
    {
    public:
        Constraint();
        ~Constraint();

        void setup(ALM *alm);

        int constraint_mode;
        int number_of_constraints;
        std::string fc2_file, fc3_file;
        bool fix_harmonic, fix_cubic;
        int constraint_algebraic;

        double **const_mat;
        double *const_rhs;
        double tolerance_constraint;

        bool exist_constraint;
        bool extra_constraint_from_symmetry;
        std::string rotation_axis;

        std::vector<ConstraintClass> *const_symmetry;
        std::vector<ConstraintTypeFix> *const_fix;
        std::vector<ConstraintTypeRelate> *const_relate;
        std::vector<ConstraintTypeRelate> *const_relate_rotation;
        boost::bimap<int, int> *index_bimap;

        void get_mapping_constraint(int,
                                    std::vector<int> *,
                                    std::vector<ConstraintClass> *,
                                    std::vector<ConstraintTypeFix> *,
                                    std::vector<ConstraintTypeRelate> *,
                                    boost::bimap<int, int> *);

    private:

        bool impose_inv_T, impose_inv_R, exclude_last_R;

        std::vector<ConstraintClass> *const_translation;
        std::vector<ConstraintClass> *const_rotation_self;
        std::vector<ConstraintClass> *const_rotation_cross;
        std::vector<ConstraintClass> *const_self;

        void set_default_variables();
        void deallocate_variables();

        int levi_civita(int,
                        int,
                        int);

        void generate_rotational_constraint(System *,
                                            Symmetry *,
                                            Interaction *,
                                            Fcs *,
                                            const int,
                                            std::vector<ConstraintClass> *,
                                            std::vector<ConstraintClass> *);

        void calc_constraint_matrix(int,
                                    std::vector<int> *,
                                    int,
                                    int &);

        void setup_rotation_axis(bool [3][3]);
        bool is_allzero(int,
                        const double *,
                        int nshift = 0);
        bool is_allzero(const std::vector<int> &,
                        int &);
        bool is_allzero(const std::vector<double> &,
                        double,
                        int &);

        void remove_redundant_rows(int,
                                   std::vector<ConstraintClass> &,
                                   double tolerance = eps12);


        void generate_symmetry_constraint_in_cartesian(int,
                                                       Symmetry *,
                                                       Interaction *,
                                                       Fcs *,
                                                       const int,
                                                       std::vector<ConstraintClass> *);

        void get_constraint_translation(const Cell &,
                                        Symmetry *,
                                        Interaction *,
                                        Fcs *,
                                        int,
                                        const std::vector<FcProperty> &,
                                        int,
                                        std::vector<ConstraintClass> &);

        void get_constraint_translation2(const Cell &,
                                        Symmetry *,
                                        Interaction *,
                                        Fcs *,
                                        int,
                                        const std::vector<FcProperty> &,
                                        int,
                                        std::vector<ConstraintClass> &);

        void generate_translational_constraint(const Cell &,
                                               Symmetry *,
                                               Interaction *,
                                               Fcs *,
                                               const int,
                                               std::vector<ConstraintClass> *);

        void fix_forceconstants_to_file(int,
                                        Symmetry *,
                                        Fcs *,
                                        std::string,
                                        std::vector<ConstraintTypeFix> &);
    };

    extern "C" {
    void dgetrf_(int *m,
                 int *n,
                 double *a,
                 int *lda,
                 int *ipiv,
                 int *info);
    }
}
