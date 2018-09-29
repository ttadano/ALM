/*
 files.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "files.h"

using namespace ALM_NS;

Files::Files()
{
    print_hessian = false;
}

Files::~Files() {}

void Files::init()
{
    file_fcs = job_title + ".fcs";
    file_hes = job_title + ".hessian";
}
