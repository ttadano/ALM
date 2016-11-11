/*
 alm.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <string>
#include "alm.h"
#include "alm_core.h"
#include "constraint.h"
#include "interaction.h"
#include "fcs.h"
#include "fitting.h"
#include "patterndisp.h"

using namespace ALM_NS;

ALM::ALM()
{
    alm_core = new ALMCore();
    alm_core->create();
}

ALM::~ALM()
{
    delete alm_core;
}

void ALM::initialize()
{
    alm_core->initialize();
}

void ALM::finalize()
{
    alm_core->finalize();
}

ALMCore * ALM::get_alm_core()
{
    return alm_core;
}

int ALM::get_fc_length(const int fc_order)  // harmonic=2, ...
{
    int id, order, num_unique_elems, num_equiv_elems;
    Fcs *fcs;

    fcs = alm_core->fcs;
    id = 0;
    order = fc_order - 2;
    if (fcs->ndup[order].size() < 1) {return 0;}
    num_unique_elems = fcs->ndup[order].size();
    for (int iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
        num_equiv_elems = fcs->ndup[order][iuniq];
        id += num_equiv_elems;
    }
    return id;
}

void ALM::get_fc(double *fc_value,
                 int *elem_indices, // (len(fc_value), fc_order) is flatten.
                 const int fc_order) // harmonic=2, ...
{
    int j, k, ip, id;
    int order, num_unique_elems, num_equiv_elems, maxorder;
    double fc_elem, coef;
    Fcs *fcs;
    Fitting *fitting;

    fcs = alm_core->fcs;
    fitting = alm_core->fitting;
    maxorder = alm_core->interaction->maxorder;
    ip = 0;
    for (order = 0; order < maxorder; ++order) {
        if (fcs->ndup[order].size() < 1) {continue;} // is this needed?
        id = 0;
        num_unique_elems = fcs->ndup[order].size();
        for (int iuniq = 0; iuniq < num_unique_elems; ++iuniq) {
            num_equiv_elems = fcs->ndup[order][iuniq];
            fc_elem = fitting->params[ip];
            for (j = 0; j < num_equiv_elems; ++j) {
                if (order == fc_order - 2) {
                    // coef is normally 1 or -1.
                    coef = fcs->fc_set[order][id].coef;
                    fc_value[id] = fc_elem * coef;
                    for (k = 0; k < fc_order; ++k) {
                        elem_indices[id * fc_order + k] = 
                            fcs->fc_set[order][id].elems[k];
                    }
                }
                ++id;
            }
            ++ip;
        }
        if (order >= fc_order - 2) {break;}
    }
}

void ALM::run_fitting()
{
    alm_core->constraint->setup();
    alm_core->fitting->fitmain();
}

void ALM::run_suggest()
{
    alm_core->displace->gen_displacement_pattern();
}


