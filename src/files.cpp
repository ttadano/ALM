/*
 files.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "files.h"
#include "error.h"
#include "memory.h"
#include "interaction.h"
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

Files::Files()
{
    file_disp_pattern = nullptr;
    print_hessian = false;
}

Files::~Files()
{
    if (file_disp_pattern) {
        deallocate(file_disp_pattern);
    }
}

void Files::init()
{
    int i;

    file_fcs = job_title + ".fcs";
    file_hes = job_title + ".hessian";

    if (alm->mode == "suggest") {

        allocate(file_disp_pattern, interaction->maxorder);

        for (i = 0; i < interaction->maxorder; ++i) {
            if (i == 0) {
                file_disp_pattern[i] = job_title + ".pattern_HARMONIC";
            } else {
                file_disp_pattern[i] = job_title + ".pattern_ANHARM"
                    + boost::lexical_cast<std::string>(i + 2);
            }
        }
    }
}
