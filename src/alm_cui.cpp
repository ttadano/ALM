/*
 alm_cui.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include "alm_core.h"
#include "alm_cui.h"
#include "input_parser.h"
#include "input_setter.h"
#include "writer.h"
#include "constraint.h"
#include "fitting.h"
#include "version.h"
#include "patterndisp.h"

using namespace ALM_NS;

ALMCUI::ALMCUI() {}

void ALMCUI::run(int narg, char **arg)
{
    // std::cout.rdbuf( NULL );
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << " +                         Program ALM                             +" << std::endl;
    std::cout << " +                             Ver.";
    std::cout << std::setw(7) << ALAMODE_VERSION;
    std::cout << "                         +" << std::endl;
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << std::endl;

    ALMCore *alm = new ALMCore();
    alm->create();

    // Here it's tricky.
    // In alm->input can access to the public variables of alm.
    // So alm->mode is set in this part.
    InputParser *parser = new InputParser();
    parser->run(alm->input, narg, arg, alm->error, alm->memory);
    delete parser;

    Writer *writer = new Writer();
    writer->write_input_vars(alm);

    alm->initialize();

    if (alm->mode == "fitting") {

        alm->constraint->setup();
        alm->fitting->fitmain();
        writer->writeall(alm);

    } else if (alm->mode == "suggest") {

        alm->displace->gen_displacement_pattern();
        writer->write_displacement_pattern(alm);

    }

    delete writer;

    alm->finalize();
    delete alm;
}

ALMCUI::~ALMCUI()
{
}

