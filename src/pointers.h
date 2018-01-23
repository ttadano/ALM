/*
 pointers.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm_core.h"

namespace ALM_NS
{
    class Pointers
    {
    public:
        Pointers(ALMCore *ptr) :
            alm(ptr),
            system(ptr->system),
            interaction(ptr->interaction),
            fcs(ptr->fcs),
            symmetry(ptr->symmetry),
            fitting(ptr->fitting),
            constraint(ptr->constraint),
            files(ptr->files),
            displace(ptr->displace),
            error(ptr->error)
        {
        }

        virtual ~Pointers()
        {
        }

    protected:
        ALMCore *alm;
        System *&system;
        Interaction *&interaction;
        Fcs *&fcs;
        Symmetry *&symmetry;
        Fitting *&fitting;
        Constraint *&constraint;
        Files *&files;
        Displace *&displace;
        Error *&error;
    };
}
