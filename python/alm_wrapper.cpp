#include "../src/alm.h"
#include "../src/memory.h"
#include "alm_wrapper.h"
#include <cstdlib>
#include <string>
#include <iostream>

extern "C" {

    #define MAX_NUM_ALM 10

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

    static ALM_NS::ALM *alm[MAX_NUM_ALM];
    static int is_alm_used[MAX_NUM_ALM];
    static int is_alm_init = 0;

    void alm_init(void)
    {
        int id;
        for (id = 0; id < MAX_NUM_ALM; id++) {
            if (alm[id]) {
                delete alm[id];
            }
            alm[id] = NULL;
            is_alm_used[id] = 0;
        }
        is_alm_init = 1;
    }

    int alm_new(void)
    {
        int id;

        if (!is_alm_init) {
            alm_init();
        }

        for (id = 0; id < MAX_NUM_ALM; id++) {
            if (!is_alm_used[id]) {
                alm[id] = new ALM_NS::ALM();
                is_alm_used[id] = 1;
                break;
            }
        }

        if (id == MAX_NUM_ALM) {
            return -1;
        } else {
            return id;
        }
    }

    void alm_delete(const int id)
    {
        if (alm[id]) {
            delete alm[id];
            alm[id] = NULL;
            is_alm_used[id] = 0;
        }
    }

    // void set_output_filename_prefix(const std::string prefix);
    // void set_is_print_symmetry(const int is_printsymmetry);
    // void set_is_print_hessians(const bool print_hessian);
    // void set_symmetry_param(const int nsym);
    // void set_symmetry_tolerance(const double tolerance);
    // void set_displacement_param(const bool trim_dispsign_for_evenfunc);
    // void set_displacement_basis(const std::string str_disp_basis);
    // void set_periodicity(const int is_periodic[3]);

    // kind_in contains integer numbers to distinguish chemical
    // elements. This is transformed to kind in ALM format, which
    // contains incrementing integer number starting from 1.
    // Here the mapping from the numbers in kind_in to those in kind
    // is made by finding unique numbers (i.e., kind_uniqe) in kind_in
    // and keeping the order, e.g., [8, 8, 4, 4] --> [1, 1, 2, 2].
    void alm_set_cell(const int id,
                      const size_t nat,
                      const double lavec[3][3],
                      const double xcoord[][3],
                      const int atomic_numbers[],
                      int kind[])
    {
        size_t i, j, nkd;
        int kind_unique[nat];
        bool in_kind_unique;

        kind_unique[0] = atomic_numbers[0];
        kind[0] = 1;
        nkd = 1;

        for (i = 1; i < nat; ++i) {
            in_kind_unique = false;
            for (j = 0; j < nkd; ++j) {
                if (kind_unique[j] == atomic_numbers[i]) {
                    in_kind_unique = true;
                    kind[i] = j + 1;
                    break;
                }
            }
            if (!in_kind_unique) {
                kind_unique[nkd] = atomic_numbers[i];
                kind[i] = nkd + 1;
                ++nkd;
            }
        }
        std::string *kdname = new std::string[nkd];
        //std::string kdname[nkd];
        for (int i = 0; i < nkd; i++) {
            kdname[i] = atom_name[abs(kind_unique[i]) % 118];
        }

        alm[id]->set_cell(nat, lavec, xcoord, kind, kdname);
        delete [] kdname;
    }

    void alm_set_verbosity(const int id, const int verbosity)
    {
        alm[id]->set_verbosity(verbosity);
    }

    // void set_magnetic_params(const unsigned int nat,
    //                          const double* const * magmom,
    //                          const bool lspin,
    //                          const int noncollinear,
    //                          const int trev_sym_mag,
    //                          const std::string str_magmom);

    void alm_set_displacement_and_force(const int id,
                                        const double* u_in,
                                        const double* f_in,
                                        const size_t nat,
                                        const size_t ndata_used)
    {
        std::vector<std::vector<double>> u, f;

        u.resize(ndata_used, std::vector<double>(3 * nat));
        f.resize(ndata_used, std::vector<double>(3 * nat));

        for (auto i = 0; i < ndata_used; i++) {
            for (auto j = 0; j < 3 * nat; j++) {
                u[i][j] = u_in[i * nat * 3 + j];
                f[i][j] = f_in[i * nat * 3 + j];
            }
        }

        alm[id]->set_training_data(u, f);

        u.clear();
        f.clear();
    }

    size_t alm_get_nrows_sensing_matrix(const int id)
    {
        return alm[id]->get_nrows_sensing_matrix();
    }

    void alm_set_constraint_type(const int id,
                                 const int constraint_flag) // ICONST
    {
        alm[id]->set_constraint_mode(constraint_flag);
    }

    // void set_fitting_constraint_rotation_axis(const std::string rotation_axis) // ROTAXIS

    void alm_define(const int id,
                    const int maxorder,
                    const size_t nkd,
                    const int *nbody_include,
                    const double *cutoff_radii_in)
    {
    /*    double ***cutoff_radii;
        int count;

        if (nkd > 0) {
            ALM_NS::allocate(cutoff_radii, maxorder, nkd, nkd);
            count = 0;
            for (size_t i = 0; i < maxorder; i++) {
                for (size_t j = 0; j < nkd; j++) {
                    for (size_t k = 0; k < nkd; k++) {
                        cutoff_radii[i][j][k] = cutoff_radii_in[count];
                        count++;
                    }
                }
            }
        } else {
            cutoff_radii = nullptr;
        }*/

        alm[id]->define(maxorder,
                        nkd,
                        nbody_include,
                        cutoff_radii_in);

 /*       if (nkd > 0) {
            ALM_NS::deallocate(cutoff_radii);
        }*/
    }

    void alm_generate_force_constant(const int id)
    {
        alm[id]->generate_force_constant();
    }

    int alm_get_atom_mapping_by_pure_translations(const int id,
                                                  int *map_p2s)
    {
        const auto map_p2s_vv = alm[id]->get_atom_mapping_by_pure_translations();

        auto nat_prim = map_p2s_vv.size();
        auto ntran = map_p2s_vv[0].size();

        size_t count = 0;

        for (size_t i = 0; i < ntran; i++) {
            for (size_t j = 0; j < nat_prim; j++) {
                map_p2s[count] = map_p2s_vv[j][i];
                count++;
            }
        }

        return ntran;
    }

    size_t alm_get_number_of_displacement_patterns(const int id,
                                                   const int fc_order) // harmonic=1,
    {
        return alm[id]->get_number_of_displacement_patterns(fc_order);
    }

    void alm_get_number_of_displaced_atoms(const int id,
                                           int *numbers,
                                           const int fc_order) // harmonic=1,
    {
        alm[id]->get_number_of_displaced_atoms(numbers, fc_order);
    }

    int alm_get_displacement_patterns(const int id,
                                      int *atom_indices,
                                      double *disp_patterns,
                                      const int fc_order) // harmonic=1,
    {
        return alm[id]->get_displacement_patterns(atom_indices,
                                                  disp_patterns,
                                                  fc_order);
    }

    size_t alm_get_number_of_fc_elements(const int id,
                                         const int fc_order)  // harmonic=1, ...
    {
        return alm[id]->get_number_of_fc_elements(fc_order);
    }

    size_t alm_get_number_of_irred_fc_elements(const int id,
                                               const int fc_order)
    {
        return alm[id]->get_number_of_irred_fc_elements(fc_order);
    }

    void alm_get_fc_origin(const int id,
                           double *fc_values,
                           int *elem_indices, // (len(fc_values), fc_order + 1) is flatten.
                           const int fc_order)
    {
        alm[id]->get_fc_origin(fc_values, elem_indices, fc_order);
    }

    void alm_get_fc_irreducible(const int id,
                                double *fc_values,
                                int *elem_indices, // (len(fc_values), fc_order + 1) is flatten.
                                const int fc_order)
    {
        alm[id]->get_fc_irreducible(fc_values, elem_indices, fc_order);
    }

    void alm_get_fc_all(const int id,
                        double *fc_values,
                        int *elem_indices, // (len(fc_values), fc_order + 1) is flatten.
                        const int fc_order)
    {
        alm[id]->get_fc_all(fc_values, elem_indices, fc_order);
    }

    void alm_set_fc(const int id, double *fc_in)
    {
        alm[id]->set_fc(fc_in);
    }

    void alm_get_matrix_elements(const int id,
                                 double *amat,
                                 double *bvec)
    {
        alm[id]->get_matrix_elements(amat, bvec);
    }

    void alm_run_suggest(const int id)
    {
        alm[id]->set_run_mode("suggest");
        alm[id]->run_suggest();
    }

    int alm_optimize(const int id, const char *solver)
    {
        alm[id]->set_run_mode("optimize");
        std::string str_solver = std::string(solver);

        int info;

        if (str_solver == "dense") {

            alm[id]->set_sparse_mode(0);
            info = alm[id]->run_optimize();

        } else if (str_solver == "SimplicialLDLT") {

            alm[id]->set_sparse_mode(1);
            info = alm[id]->run_optimize();

        } else {
            std::cerr << " Unsupported solver type : " << str_solver << std::endl;
            return EXIT_FAILURE;
        }
        return info;
    }
}
