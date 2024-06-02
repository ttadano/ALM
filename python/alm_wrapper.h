#ifndef __ALM_WRAPPER_H__
#define __ALM_WRAPPER_H__

#ifdef __cplusplus
extern "C" {
#endif

struct optimizer_control {
    int linear_model;      // 1 : least-squares, 2 : elastic net
    int use_sparse_solver; // 0: No, 1: Yes
    int maxnum_iteration;
    double tolerance_iteration;
    int output_frequency;

    // Options related to L1-regularized optimization
    int standardize;
    double displacement_normalization_factor;
    int debiase_after_l1opt;

    // cross-validation related variables
    int cross_validation; // 0 : No CV mode, -1 or > 0: CV mode
    double l1_alpha;      // L1-regularization coefficient
    double l1_alpha_min;
    double l1_alpha_max;
    int num_l1_alpha;
    double l1_ratio; // l1_ratio = 1 for LASSO; 0 < l1_ratio < 1 for Elastic net
    int save_solution_path;
    int mirror_image_conv;
};

void alm_init(void);
int alm_new(void);
void alm_delete(const int id);
void alm_define(const int id,
                const int maxorder,
                const size_t nkd,
                const int *nbody_include,
                const double *cutoff_radii_in,
                const char *fc_basis);
void alm_suggest(const int id);
int alm_optimize(const int id,
                 const char *solver);
void alm_init_fc_table(const int id);
void alm_set_optimizer_control(const int id,
                               const struct optimizer_control optcontrol,
                               const int updated[15]);
void alm_set_cell(const int id,
                  const size_t nat,
                  const double lavec[3][3],
                  const double xcoord[][3],
                  const int numbers[],
                  const size_t nkind,
                  const int kind_numbers[]);
void alm_set_verbosity(const int id,
                       const int verbosity);
void alm_set_u_train(const int id,
                     const double *u_in,
                     const size_t nat,
                     const size_t ndata_used);
void alm_set_f_train(const int id,
                     const double *f_in,
                     const size_t nat,
                     const size_t ndata_used);
void alm_set_constraint_type(const int id,
                             const int constraint_flag); // ICONST
void alm_set_forceconstants_to_fix(const int id,
                                   const int *fc_indices,
                                   const double *fc_values,
                                   const size_t nfcs,
                                   const int fc_order); // harmonic=1
void alm_set_fc(const int id, double *fc_in);

void alm_set_output_filename_prefix(const int id,
                                    const char *prefix_in);
// void set_magnetic_params(const double* const * magmom,
//                         const bool lspin,
//                         const int noncollinear,
//                         const int trev_sym_mag,
//                         const std::string str_magmom);
// void set_fitting_constraint_rotation_axis(const std::string rotation_axis) // ROTAXIS
// void set_fitting_filenames(const std::string dfile,
//                           const std::string ffile);
struct optimizer_control alm_get_optimizer_control(const int id);
int alm_get_u_train(const int id, double *u_out, const int nelems_in);
int alm_get_f_train(const int id, double *f_out, const int nelems_in);
double alm_get_cv_l1_alpha(const int id);
int alm_get_atom_mapping_by_pure_translations(const int id,
                                              int *map_p2s);
size_t alm_get_number_of_displacement_patterns(const int id,
                                               const int fc_order); // harmonic=1,
void alm_get_number_of_displaced_atoms(const int id,
                                       int *numbers,
                                       const int fc_order); // harmonic=1,
size_t alm_get_number_of_data(const int id);
size_t alm_get_nrows_sensing_matrix(const int id);
int alm_get_displacement_patterns(const int id,
                                  int *atom_indices,
                                  double *disp_patterns,
                                  const int fc_order); // harmonic=1,
size_t alm_get_number_of_fc_elements(const int id,
                                     const int fc_order); // harmonic=1, ...
size_t alm_get_number_of_fc_origin(const int id,
                                   const int fc_order, // harmonic=1, ...
                                   const int permutation);
size_t alm_get_number_of_irred_fc_elements(const int id,
                                           const int fc_order); // harmonic=1, ...
void alm_get_fc_origin(const int id,
                       double *fc_value,
                       int *elem_indices, // (len(fc_value), fc_order + 1) is flatten.
                       const int fc_order,
                       const int permutation);
void alm_get_fc_irreducible(const int id,
                            double *fc_value,
                            int *elem_indices, // (len(fc_value), fc_order + 1) is flatten.
                            const int fc_order);
void alm_get_fc_all(const int id,
                    double *fc_value,
                    int *elem_indices, // (len(fc_value), fc_order + 1) is flatten.
                    const int fc_order,
                    const int permutation);
void alm_get_matrix_elements(const int id,
                             double *amat,
                             double *bvec);

void alm_save_fc(const int id,
                 const char *filename,
                 const char *format,
                 const int maxorder_to_save);

void alm_get_fc_dependency(const int id,
                           int *elem_indices_irred,
                           int *elem_indices_orig,
                           double *dependency_mat,
                           const int fc_order);


#ifdef __cplusplus
}
#endif

#endif
