/*
 alm_cui.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include "alm.h"
#include "alm_core.h"
#include "alm_cui.h"
#include "input_parser.h"
#include "writer.h"
#include "version.h"
#include "timer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace ALM_NS;

ALMCUI::ALMCUI() {}

void ALMCUI::run(int narg, char **arg)
{
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << " +                         Program ALM                             +" << std::endl;
    std::cout << " +                             Ver.";
    std::cout << std::setw(7) << ALAMODE_VERSION;
    std::cout << "                         +" << std::endl;
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << std::endl;

    ALM *alm = new ALM();

#ifdef _OPENMP
    std::cout << " Number of OpenMP threads = "
        << omp_get_max_threads() << std::endl << std::endl;
#endif

    ALMCore *alm_core = alm->get_alm_core();
    // alm_core->mode is set herein.
    InputParser *input_parser = new InputParser();
    input_parser->run(alm_core, narg, arg);

    Writer *writer = new Writer();
    writer->write_input_vars(alm);


    if (alm_core->mode == "fitting") {
        input_parser->parse_displacement_and_force(alm_core);
    }
    delete input_parser;

    alm->run();

    if (alm_core->mode == "fitting") {
        writer->writeall(alm);
    } else if (alm_core->mode == "suggest") {
        writer->write_displacement_pattern(alm);
    }
    delete writer;

    delete alm;
}

ALMCUI::~ALMCUI()
{
}

