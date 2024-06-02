/*
 alm_cui.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "alm_cui.h"
#include "alm.h"
#include "input_parser.h"
#include "timer.h"
#include "version.h"
#include <iostream>
#include <iomanip>

#ifdef _OPENMP

#include <omp.h>

#endif

using namespace ALM_NS;

ALMCUI::ALMCUI() {}

ALMCUI::~ALMCUI() {}

void ALMCUI::run(const int narg,
                 char **arg) const
{
    auto alm = new ALM();

    // alm->mode is set herein.
    auto input_parser = new InputParser();
    input_parser->run(alm, narg, arg);
    auto run_mode = input_parser->get_run_mode();
    delete input_parser;

    if (alm->get_verbosity() > 0) {
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << " +                         Program ALM                             +" << std::endl;
        std::cout << " +                             Ver.";
        std::cout << std::setw(7) << ALAMODE_VERSION;
        std::cout << "                         +" << std::endl;
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << std::endl;
#ifdef _OPENMP
        std::cout << " Number of OpenMP threads = "
                  << omp_get_max_threads() << std::endl << std::endl;
#endif

        std::cout << " Job started at " << alm->timer->DateAndTime() << std::endl;
    }

    if (alm->get_verbosity() > 0) {
        alm->writer->write_input_vars(alm->system,
                                      alm->symmetry,
                                      alm->cluster,
                                      alm->displace,
                                      alm->fcs,
                                      alm->constraint,
                                      alm->optimize,
                                      alm->files,
                                      run_mode);
    }

    alm->init_fc_table();

    if (run_mode == "optimize") {
        alm->run_optimize();
        if (alm->get_optimizer_control().linear_model == 1 ||
            (alm->get_optimizer_control().linear_model >= 2
             && alm->get_optimizer_control().cross_validation == 0)) {
            alm->writer->writeall(alm->system,
                                  alm->symmetry,
                                  alm->cluster,
                                  alm->constraint,
                                  alm->fcs,
                                  alm->optimize,
                                  alm->files,
                                  alm->get_verbosity());

            const auto nfcs_irred = alm->get_number_of_irred_fc_elements(2);
            const auto nfcs_origin = alm->get_number_of_fc_origin(2, 0);

            std::cout << nfcs_irred << " " << nfcs_origin << '\n';

            double *flattened_array = new double[nfcs_irred * nfcs_origin];
            int *index_elements_origin = new int[nfcs_origin * 3];
            int *index_elements_irred = new int[nfcs_irred * 3];

            alm->get_fc_dependency_mat(2,
                                       index_elements_irred,
                                       index_elements_origin,
                                       flattened_array);

            for (size_t i = 0; i < nfcs_origin; ++i) {
                for (size_t j = 0; j < nfcs_irred; ++j) {
                    std::cout << flattened_array[i * nfcs_irred + j] << " ";
                }
                std::cout << '\n';
            }

            delete[] flattened_array;
            delete[] index_elements_irred;
            delete[] index_elements_origin;

        }
    } else if (run_mode == "suggest") {
        alm->run_suggest();
        alm->writer->write_displacement_pattern(alm->cluster,
                                                alm->displace,
                                                alm->files->get_prefix(),
                                                alm->get_verbosity());
    }

    if (alm->get_verbosity() > 0) {
        std::cout << std::endl << " Job finished at "
                  << alm->timer->DateAndTime() << std::endl;
    }

    delete alm;
}
