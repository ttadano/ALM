#ifndef __ALM_WRAPPER_H__
#define __ALM_WRAPPER_H__

#ifdef __cplusplus
extern "C" {
#endif

    void alm_init(void);
    int alm_new(void);
    void alm_delete(const int id);
    // void set_output_filename_prefix(const std::string prefix);
    // void set_is_print_symmetry(const int is_printsymmetry);
    // void set_is_print_hessians(const bool print_hessian);
    // void set_symmetry_param(const int nsym);
    // void set_symmetry_tolerance(const double tolerance);
    // void set_displacement_param(const bool trim_dispsign_for_evenfunc);
    // void set_displacement_basis(const std::string str_disp_basis);
    // void set_periodicity(const int is_periodic[3]);
    void alm_set_cell(const int id,
                      const size_t nat,
                      const double lavec[3][3],
                      const double xcoord[][3],
                      const int kd[],
                      int kind[]);
    void alm_set_verbosity(const int id,
                           const int verbosity);
    // void set_magnetic_params(const double* const * magmom,
    //                         const bool lspin,
    //                         const int noncollinear,
    //                         const int trev_sym_mag,
    //                         const std::string str_magmom);
    void alm_set_displacement_and_force(const int id,
                                        const double* u_in,
                                        const double* f_in,
                                        const size_t nat,
                                        const size_t ndata_used);

    void alm_set_constraint_type(const int id,
                                 const int constraint_flag); // ICONST
    // void set_fitting_constraint_rotation_axis(const std::string rotation_axis) // ROTAXIS
    // void set_fitting_filenames(const std::string dfile,
    //                           const std::string ffile);
    void alm_define(const int id,
                    const int maxorder,
                    const size_t nkd,
                    const int *nbody_include,
                    const double *cutoff_radii);
    void alm_generate_force_constant(const int id);
    int alm_get_atom_mapping_by_pure_translations(const int id,
                                                  int *map_p2s);
    size_t alm_get_number_of_displacement_patterns(const int id,
                                                   const int fc_order); // harmonic=1,
    void alm_get_number_of_displaced_atoms(const int id,
                                           int *numbers,
                                           const int fc_order); // harmonic=1,
    int alm_get_displacement_patterns(const int id,
                                      int *atom_indices,
                                      double *disp_patterns,
                                      const int fc_order); // harmonic=1,
    size_t alm_get_number_of_fc_elements(const int id,
                                         const int fc_order); // harmonic=1, ...

    size_t alm_get_number_of_irred_fc_elements(const int id,
                                               const int fc_order); // harmonic=1, ...

    void alm_get_fc_origin(const int id,
                    double *fc_value,
                    int *elem_indices, // (len(fc_value), fc_order + 1) is flatten.
                    const int fc_order);

    void alm_get_fc_irreducible(const int id,
                    double *fc_value,
                    int *elem_indices, // (len(fc_value), fc_order + 1) is flatten.
                    const int fc_order);

    void alm_get_fc_all(const int id,
                    double *fc_value,
                    int *elem_indices, // (len(fc_value), fc_order + 1) is flatten.
                    const int fc_order);

    void alm_set_fc(const int id, double *fc_in);

    void alm_get_matrix_elements(const int id,
                                 double *amat,
                                 double *bvec);
    size_t alm_get_nrows_sensing_matrix(const int id);

    void alm_run_suggest(const int id);
    int alm_optimize(const int id,
                     const char *solver);

#ifdef __cplusplus
}
#endif

#endif
