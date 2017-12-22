/*
 alm_core.h

 Copyright (c) 2014 2015 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

/* Declaration of pointers used in the whole program. */

#pragma once
#include <string>

namespace ALM_NS
{
    class ALMCore
    {
    public:
        class InputSetter *input;
        class Memory *memory;
        class System *system;
        class Interaction *interaction;
        class Fcs *fcs;
        class Symmetry *symmetry;
        class Fitting *fitting;
        class Constraint *constraint;
        class Files *files;
        class Displace *displace;
        class Error *error;
        class Timer *timer;
        ALMCore();
        ~ALMCore();
        void create();
        void initialize();
        void finalize();

        std::string mode;
        bool print_hessian;
    };
}
