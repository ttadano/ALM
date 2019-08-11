#include "../src/alm.h"
#include "../src/memory.h"
#include "../src/optimize.h"
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

    void alm_define(const int id,
                    const int maxorder,
                    const size_t nkd,
                    const int *nbody_include,
                    const double *cutoff_radii_in,
                    const char *fc_basis)
    {
        std::string str_fc_basis = std::string(fc_basis);

        if (str_fc_basis == "Lattice" || str_fc_basis == "Cartesian") {
            alm[id]->set_forceconstant_basis(str_fc_basis);
        }

        alm[id]->define(maxorder,
                        nkd,
                        nbody_include,
                        cutoff_radii_in);
    }

    void alm_suggest(const int id)
    {
        alm[id]->run_suggest();
    }

    int alm_optimize(const int id, const char *solver)
    {
        std::string str_solver = std::string(solver);

        int info;

        auto optctrl = alm[id]->get_optimizer_control();
        if (str_solver == "dense") {
            optctrl.use_sparse_solver = 0;
        } else if (str_solver == "SimplicialLDLT") {
            optctrl.use_sparse_solver = 1;
        } else {
            std::cerr << " Unsupported solver type : " << str_solver << std::endl;
            return EXIT_FAILURE;
        }
        alm[id]->set_optimizer_control(optctrl);

        info = alm[id]->run_optimize();

        return info;
    }

    void alm_init_fc_table(const int id)
    {
        alm[id]->init_fc_table();
    }

    void alm_set_optimizer_control(const int id,
                                   const struct optimizer_control optcontrol,
                                   const int updated[15])
    {
        auto optctrl = alm[id]->get_optimizer_control();

        if (updated[0]) {
            optctrl.linear_model = optcontrol.linear_model;
        }
        if (updated[1]) {
            optctrl.use_sparse_solver = optcontrol.use_sparse_solver;
        }
        if (updated[2]) {
            optctrl.maxnum_iteration = optcontrol.maxnum_iteration;
        }
        if (updated[3]) {
            optctrl.tolerance_iteration = optcontrol.tolerance_iteration;
        }
        if (updated[4]) {
            optctrl.output_frequency = optcontrol.output_frequency;
        }
        if (updated[5]) {
            optctrl.standardize = optcontrol.standardize;
        }
        if (updated[6]) {
            optctrl.displacement_normalization_factor = optcontrol.displacement_normalization_factor;
        }
        if (updated[7]) {
            optctrl.debiase_after_l1opt = optcontrol.debiase_after_l1opt;
        }
        if (updated[8]) {
            optctrl.cross_validation = optcontrol.cross_validation;
        }
        if (updated[9]) {
            optctrl.l1_alpha = optcontrol.l1_alpha;
        }
        if (updated[10]) {
            optctrl.l1_alpha_min = optcontrol.l1_alpha_min;
        }
        if (updated[11]) {
            optctrl.l1_alpha_max = optcontrol.l1_alpha_max;
        }
        if (updated[12]) {
            optctrl.num_l1_alpha = optcontrol.num_l1_alpha;
        }
        if (updated[13]) {
            optctrl.l1_ratio = optcontrol.l1_ratio;
        }
        if (updated[14]) {
            optctrl.save_solution_path = optcontrol.save_solution_path;
        }
        alm[id]->set_optimizer_control(optctrl);
    }

    void alm_set_cell(const int id,
                      const size_t nat,
                      const double lavec[3][3],
                      const double xcoord[][3],
                      const int numbers[],
                      const size_t nkind,
                      const int kind_numbers[])
    {
        int atom_mapping[nat];
        std::string *kdname = new std::string[nkind];

        for (auto i = 0; i < nkind; i++) {
            kdname[i] = atom_name[kind_numbers[i]];
        }

        for (auto i = 0; i < nat; i++) {
            for (auto j = 0; j < nkind; j++) {
                if (numbers[i] == kind_numbers[j]) {
                    atom_mapping[i] = j + 1;
                    break;
                }
            }
        }

        alm[id]->set_cell(nat, lavec, xcoord, atom_mapping, kdname);
        delete [] kdname;
    }

    void alm_set_verbosity(const int id, const int verbosity)
    {
        alm[id]->set_verbosity(verbosity);
    }

    // numbers: atomic numbers
    // kind_numbers: unique atomic numbers preserving appering order
    // atom_mapping: numbers mapped from atomic numbers. Each number
    //               corresponds to the index of  the atomic number
    //               found in kind_numbers. Note this index starts
    //               with 1.
    // kdname: List of element names corresponds to that of kind_numbers.

    // void set_magnetic_params(const unsigned int nat,
    //                          const double* const * magmom,
    //                          const bool lspin,
    //                          const int noncollinear,
    //                          const int trev_sym_mag,
    //                          const std::string str_magmom);

    void alm_set_u_train(const int id,
                         const double* u_in,
                         const size_t nat,
                         const size_t ndata_used)
    {
        std::vector<std::vector<double>> u;

        u.resize(ndata_used, std::vector<double>(3 * nat));

        for (size_t i = 0; i < ndata_used; i++) {
            for (size_t j = 0; j < 3 * nat; j++) {
                u[i][j] = u_in[i * nat * 3 + j];
            }
        }

        alm[id]->set_u_train(u);

        u.clear();
    }

    void alm_set_f_train(const int id,
                         const double* f_in,
                         const size_t nat,
                         const size_t ndata_used)
    {
        std::vector<std::vector<double>> f;

        f.resize(ndata_used, std::vector<double>(3 * nat));

        for (size_t i = 0; i < ndata_used; i++) {
            for (size_t j = 0; j < 3 * nat; j++) {
                f[i][j] = f_in[i * nat * 3 + j];
            }
        }

        alm[id]->set_f_train(f);

        f.clear();
    }

    void alm_set_constraint_type(const int id,
                                 const int constraint_flag) // ICONST
    {
        alm[id]->set_constraint_mode(constraint_flag);
    }

    void alm_set_fc(const int id, double *fc_in)
    {
        alm[id]->set_fc(fc_in);
    }

    void alm_set_output_filename_prefix(const int id,
                                        const char *prefix_in) {
        std::string prefix(prefix_in);
        alm[id]->set_output_filename_prefix(prefix);
    }

    struct optimizer_control alm_get_optimizer_control(const int id)
    {
        struct optimizer_control optcontrol;
        auto optctrl = alm[id]->get_optimizer_control();

        optcontrol.linear_model = optctrl.linear_model;
        optcontrol.use_sparse_solver = optctrl.use_sparse_solver;
        optcontrol.maxnum_iteration = optctrl.maxnum_iteration;
        optcontrol.tolerance_iteration = optctrl.tolerance_iteration;
        optcontrol.output_frequency = optctrl.output_frequency;
        optcontrol.standardize = optctrl.standardize;
        optcontrol.displacement_normalization_factor = optctrl.displacement_normalization_factor;
        optcontrol.debiase_after_l1opt = optctrl.debiase_after_l1opt;
        optcontrol.cross_validation = optctrl.cross_validation;
        optcontrol.l1_alpha = optctrl.l1_alpha;
        optcontrol.l1_alpha_min = optctrl.l1_alpha_min;
        optcontrol.l1_alpha_max = optctrl.l1_alpha_max;
        optcontrol.num_l1_alpha = optctrl.num_l1_alpha;
        optcontrol.l1_ratio = optctrl.l1_ratio;
        optcontrol.save_solution_path = optctrl.save_solution_path;

        return optcontrol;
    }

    int alm_get_u_train(const int id, double* u_out, const int nelems_in)
    {
        auto u = alm[id]->get_u_train();

        // Data is not copied.
        if (u.size() * u[0].size() != nelems_in) {
            return 0;
        }

        for (size_t i = 0; i < u.size(); i++) {
            for (size_t j = 0; j < u[0].size(); j++) {
                u_out[i * u[0].size() + j] = u[i][j];
            }
        }

        // Succeeded
        return 1;
    }

    int alm_get_f_train(const int id, double* f_out, const int nelems_in)
    {
        auto f = alm[id]->get_f_train();

        // Data is not copied.
        if (f.size() * f[0].size() != nelems_in) {
            return 0;
        }

        for (size_t i = 0; i < f.size(); i++) {
            for (size_t j = 0; j < f[0].size(); j++) {
                f_out[i * f[0].size() + j] = f[i][j];
            }
        }

        // Succeeded
        return 1;
    }

    double alm_get_cv_l1_alpha(const int id)
    {
        return alm[id]->get_cv_l1_alpha();
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

    size_t alm_get_number_of_data(const int id)
    {
        return alm[id]->get_number_of_data();
    }

    size_t alm_get_nrows_sensing_matrix(const int id)
    {
        return alm[id]->get_nrows_sensing_matrix();
    }

    // void set_fitting_constraint_rotation_axis(const std::string rotation_axis) // ROTAXIS

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

    size_t alm_get_number_of_fc_origin(const int id,
                                       const int fc_order, // harmonic = 1
                                       const int permutation)
    {
        return alm[id]->get_number_of_fc_origin(fc_order, permutation);
    }

    size_t alm_get_number_of_irred_fc_elements(const int id,
                                               const int fc_order)
    {
        return alm[id]->get_number_of_irred_fc_elements(fc_order);
    }

    void alm_get_fc_origin(const int id,
                           double *fc_values,
                           int *elem_indices, // (len(fc_values), fc_order + 1) is flatten.
                           const int fc_order,
                           const int permutation)
    {
        alm[id]->get_fc_origin(fc_values, elem_indices, fc_order, permutation);
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
                        const int fc_order,
                        const int permutation)
    {
        alm[id]->get_fc_all(fc_values, elem_indices, fc_order, permutation);
    }

    void alm_get_matrix_elements(const int id,
                                 double *amat,
                                 double *bvec)
    {
        alm[id]->get_matrix_elements(amat, bvec);
    }

}
