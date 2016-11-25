/*
 alm_core.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iomanip>
#include "alm_core.h"
#include "interaction.h"
#include "symmetry.h"
#include "system.h"
#include "files.h"
#include "fcs.h"
#include "fitting.h"
#include "constraint.h"
#include "timer.h"
#include "patterndisp.h"
#include "error.h"


using namespace ALM_NS;

ALMCore::ALMCore()
{
    // Default values
    mode = "suggest";
}

ALMCore::~ALMCore()
{
    delete files;
    delete interaction;
    delete fcs;
    delete symmetry;
    delete system;
    delete fitting;
    delete constraint;
    delete displace;
    delete error;
    delete timer;
}

void ALMCore::create()
{
    files = new Files(this);
    system = new System(this);
    interaction = new Interaction(this);
    fcs = new Fcs(this);
    symmetry = new Symmetry(this);
    fitting = new Fitting(this);
    constraint = new Constraint(this);
    displace = new Displace(this);
    error = new Error(this);
    timer = new Timer();
}

void ALMCore::initialize()
{
    system->init();
    files->init();
    symmetry->init();
    interaction->init();
    fcs->init();
}
