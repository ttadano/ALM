#pragma once

#include <vector>

void remove_redundant_rows(int,
                           std::vector<std::vector<double>> &,
                           double tolerance = 1.0e-12);

void rref(int,
          int,
          double **,
          int &,
          double tolerance = 1.0e-12);

void rref(std::vector<std::vector<double>> &,
          double tolerance = 1.0e-12);
