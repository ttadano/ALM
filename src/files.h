/*
 files.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <fstream>

namespace ALM_NS
{
    class Files
    {
    public:
        Files();
        ~Files();

        void init();

        bool print_hessian;
        std::string file_fcs, file_hes;
        std::string file_disp, file_force;

        void set_prefix(const std::string);
        std::string get_prefix() const;

    private:

        std::string job_title;
    };
}
