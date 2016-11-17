/*
 alm_core.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include "alm_core.h"
#include "interaction.h"
#include "symmetry.h"
#include "system.h"
#include "files.h"
#include "memory.h"
#include "fcs.h"
#include "fitting.h"
#include "constraint.h"
#include "timer.h"
#include "patterndisp.h"
#include "error.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace ALM_NS;

ALMCore::ALMCore()
{
#ifdef _OPENMP
    std::cout << " Number of OpenMP threads = " 
              << omp_get_max_threads() << std::endl << std::endl;
#endif

    timer = new Timer(this);

    std::cout << " Job started at " << timer->DateAndTime() << std::endl;

    // Default values
    mode = "suggest";
    print_hessian = false;
}

void ALMCore::create()
{
    memory = new Memory(this);
    files = new Files(this);
    system = new System(this);
    interaction = new Interaction(this);
    fcs = new Fcs(this);
    symmetry = new Symmetry(this);
    fitting = new Fitting(this);
    constraint = new Constraint(this);
    displace = new Displace(this);
    error = new Error(this);
}

void ALMCore::initialize()
{
    system->init();
    files->init();
    symmetry->init();
    interaction->init();
    fcs->init();
}

ALMCore::~ALMCore()
{
    std::cout << std::endl << " Job finished at " 
        << timer->DateAndTime() << std::endl;
    delete timer;
    delete files;
    delete interaction;
    delete fcs;
    delete symmetry;
    delete system;
    delete fitting;
    delete constraint;
    delete displace;
    delete memory;
    delete error;
}

