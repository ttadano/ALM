/*
 main.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <stdlib.h>
#include "alm_core.h"
#include "alm_cui.h"

using namespace ALM_NS;

int main(int argc, char **argv)
{
    ALMCore *alm = new ALMCore();
    ALMCUI *almcui = new ALMCUI();

    almcui->run(alm, argc, argv);

    delete almcui;
    delete alm;

    return EXIT_SUCCESS;
}
