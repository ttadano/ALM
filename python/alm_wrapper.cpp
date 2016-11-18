#include "alm.h"
#include "alm_wrapper.h"

extern "C" {
    static const std::string[] = {
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

    int alm_get_fc(double *fc_value,
                   int *elem_indices, // (len(fc_value), fc_order + 1) is flatten.
                   const int fc_order)
    {
    }

    int alm_get_fc_length(const int fc_order)  // harmonic=1, ...
    {
    }

    void alm_set_cell(const int nat,
                      const double lavec[3][3],
                      const double xcoord[][3],
                      const int kd[])
    {
        int i, j, nkd;
        int nkd_vals[nat];
        bool kd_exist;

        nkd_vals[0] = kd[0];
        nkd = 1;
        for (i = 1; i < nat; ++i) {
            kd_exist = false;
            for (j = 0; j < nkd; ++j) {
                if (nkd_vals[j] == kd[i]) {
                    kd_exist = true;
                    break;
                }
            }
            if (!kd_exist) {
                nkd_vals[nkd] = kd[i];
                ++nkd;
            }
        }
        std::string kdname[nkd];
        for (int i = 0; i < nkd; i++) {
            kdname[i] = abs(nkd_vals[i]) % 118;
        }
        alm->set_cell(nat, lavec, xcoord, kd, kdname);
    }

    void alm_set_interaction_vars(const int maxorder, // NORDER harmonic=1
                                  const int *nbody_include) // NBODY
    {
    }

    void alm_set_cutoff_radii(const double * rcs)
    {
    }

}
