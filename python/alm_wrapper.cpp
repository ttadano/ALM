#include "../src/alm.h"
#include "alm_wrapper.h"
#include <cstdlib>
#include <string>
#include <iostream>

extern "C" {
    static const std::string atom_name[] = {
        "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F",
        "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K",
        "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
        "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
        "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
        "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
        "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
        "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
        "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
        "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"};

    static ALM_NS::ALM *alm;

    void alm_new(void)
    {
        alm = new ALM_NS::ALM();
    }

    void alm_delete(void)
    {
        delete alm;
    }

    // void set_output_filename_prefix(const std::string prefix);
    // void set_is_print_symmetry(const int is_printsymmetry);
    // void set_is_print_hessians(const bool print_hessian);
    // void set_symmetry_param(const int nsym);
    // void set_symmetry_tolerance(const double tolerance);
    // void set_displacement_param(const bool trim_dispsign_for_evenfunc);
    // void set_displacement_basis(const std::string str_disp_basis);
    // void set_periodicity(const int is_periodic[3]);

    void alm_set_cell(const int nat,
                      const double lavec[3][3],
                      const double xcoord[][3],
                      const int kd[])
    {
        int i, j, nkd;
        int nkd_vals[nat], kd_new[nat];
        bool kd_exist;

        nkd_vals[0] = kd[0];
        kd_new[0] = 1;
        nkd = 1;
        for (i = 1; i < nat; ++i) {
            kd_exist = false;
            for (j = 0; j < nkd; ++j) {
                if (nkd_vals[j] == kd[i]) {
                    kd_exist = true;
                    kd_new[i] = j + 1;
                    break;
                }
            }
            if (!kd_exist) {
                nkd_vals[nkd] = kd[i];
                kd_new[i] = nkd + 1;
                ++nkd;
            }
        }
        std::string kdname[nkd];
        for (int i = 0; i < nkd; i++) {
            kdname[i] = atom_name[abs(nkd_vals[i]) % 118];
        }


        alm->set_cell(nat, lavec, xcoord, kd_new, kdname);
    }

    // void set_magnetic_params(const double* magmom,
    //   		       const bool lspin,
    //   		       const int noncollinear,
    //   		       const int trev_sym_mag,
    //   		       const std::string str_magmom);

    void alm_set_displacement_and_force(const double* u_in,
                                        const double* f_in,
                                        const int nat,
                                        const int ndata_used)
    {
        alm->set_displacement_and_force(u_in, f_in, nat, ndata_used);
    }

    // void set_fitting_constraint(const int constraint_flag,
    //   			  const std::string rotation_axis);
    // void set_multiplier_option(const int multiply_data);
    // void set_fitting_filenames(const std::string dfile,
    //   			 const std::string ffile);
    void alm_set_norder(const int maxorder)
    {
        alm->set_norder(maxorder);
    }

    void alm_set_interaction_range(const int *nbody_include)
    {
        alm->set_interaction_range(nbody_include);
    }

    void alm_set_cutoff_radii(const double * rcs)
    {
        alm->set_cutoff_radii(rcs);
    }

    int alm_get_number_of_displacement_patterns(const int fc_order) // harmonic=1,
    {
        return alm->get_number_of_displacement_patterns(fc_order);
    }

    void alm_get_numbers_of_displacements(int *numbers,
                                          const int fc_order) // harmonic=1,
    {
        alm->get_numbers_of_displacements(numbers, fc_order);        
    }

    int alm_get_displacement_patterns(int *atom_indices,
                                      double *disp_patterns,
                                      const int fc_order) // harmonic=1,
    {
        return alm->get_displacement_patterns(atom_indices,
                                              disp_patterns,
                                              fc_order);
    }

    int alm_get_number_of_fc_elements(const int fc_order)  // harmonic=1, ...
    {
        return alm->get_number_of_fc_elements(fc_order);
    }

    void alm_get_fc(double *fc_values,
                    int *elem_indices, // (len(fc_values), fc_order + 1) is flatten.
                    const int fc_order)
    {
        alm->get_fc(fc_values, elem_indices, fc_order);
    }

    void alm_run_suggest(void)
    {
        alm->set_run_mode("suggest");
        alm->run();
    }

    void alm_run_fitting(void)
    {
        alm->set_run_mode("fitting");
        alm->run();
    }
}
