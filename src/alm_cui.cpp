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
#include "writes.h"
#include "constraint.h"
#include "fitting.h"
#include "version.h"
#include "patterndisp.h"

using namespace ALM_NS;

ALMCUI::ALMCUI() {}

void ALMCUI::run(ALMCore *alm, int narg, char **arg)
{
    // std::cout.rdbuf( NULL );
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << " +                         Program ALM                             +" << std::endl;
    std::cout << " +                             Ver.";
    std::cout << std::setw(7) << ALAMODE_VERSION;
    std::cout << "                         +" << std::endl;
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << std::endl;

    alm->create();


    InputParser *parser = new InputParser();
    parser->parse_input(alm->input, narg, arg, alm->error, alm->memory, alm->mode);
    delete parser;

    alm->writes->write_input_vars();
    alm->initialize();

    if (alm->mode == "fitting") {

        alm->constraint->setup();
        alm->fitting->fitmain();
        alm->writes->writeall();

    } else if (alm->mode == "suggest") {

        alm->displace->gen_displacement_pattern();
        alm->writes->write_displacement_pattern();

    }

    alm->finalize();
}

ALMCUI::~ALMCUI()
{
}

