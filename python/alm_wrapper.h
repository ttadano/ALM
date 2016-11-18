#ifndef __ALM_WRAPPER_H__
#define __ALM_WRAPPER_H__

#ifdef __cplusplus 
extern "C" { 
#endif

    void alm_new(void);
    void alm_delete(void);
    // void set_output_filename_prefix(const std::string prefix);
    // void set_is_print_symmetry(const int is_printsymmetry);
    // void set_is_print_hessians(const bool print_hessian);
    // void set_symmetry_params(const int nsym,
    //   		       const double tolerance);
    // void set_displacement_params(const std::string str_disp_basis,
    //   			   const bool trim_dispsign_for_evenfunc);
    // void set_periodicity(const int is_periodic[3]);
    void alm_set_cell(const int nat,
                      const double lavec[3][3],
                      const double xcoord[][3],
                      const int kd[]);
    // void set_magnetic_params(const double* magmom,
    //   		       const bool lspin,
    //   		       const int noncollinear,
    //   		       const int trev_sym_mag,
    //   		       const std::string str_magmom);
    void alm_set_displacement_and_force(const double* u_in,
                                        const double* f_in,
                                        const int nat,
                                        const int ndata_used);
    // void set_fitting_constraint(const int constraint_flag,
    //   			  const std::string rotation_axis);
    // void set_multiplier_option(const int multiply_data);
    // void set_fitting_filenames(const std::string dfile,
    //   			 const std::string ffile);
    void alm_set_norder(const int maxorder);
    void alm_set_interaction_range(const int *nbody_include);
    void alm_set_cutoff_radii(const double * rcs);
    int alm_get_fc_length(const int fc_order);  // harmonic=1, ...
    void alm_get_fc(double *fc_value,
                    int *elem_indices, // (len(fc_value), fc_order + 1) is flatten.
                    const int fc_order);
    void alm_run_suggest(void);
    void alm_run_fitting(void);


#ifdef __cplusplus 
}
#endif

#endif
