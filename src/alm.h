/*
 alm.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm_core.h"
#include "fcs.h"
#include "fitting.h"

namespace ALM_NS
{
    class ALM
    {
    public:
        ALM();
        ~ALM();

	void initialize();
	void finalize();
	ALMCore * get_alm_core();
	int get_fc_length(const int fc_order);  // harmonic=2, ...
	void get_fc(double *fc_value,
		    int *elem_indices, // (len(fc_value), fc_order) is flatten.
		    const int fc_order); // harmonic=2, ...
	void run_fitting();
	void run_suggest();

    private:
	ALMCore *alm_core;
    };
}
