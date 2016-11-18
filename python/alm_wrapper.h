#ifndef __ALM_WRAPPER_H__
#define __ALM_WRAPPER_H__

#ifdef __cplusplus 
extern "C" { 
#endif
    void alm_new(void);
    void alm_delete(void);
    int alm_get_fc(double *fc_value,
                   int *elem_indices, // (len(fc_value), fc_order + 1) is flatten.
                   const int fc_order);
    int alm_get_fc_length(const int fc_order);  // harmonic=1, ...
    void alm_set_cell(const int nat,
                      const double lavec[3][3],
                      const double xcoord[][3],
                      const int kd[]);
    void alm_set_interaction_vars(const int maxorder, // NORDER harmonic=1
                                  const int *nbody_include); // NBODY
    void alm_set_cutoff_radii(const double * rcs);
#ifdef __cplusplus 
}
#endif

#endif
