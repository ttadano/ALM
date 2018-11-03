#pragma once

#include "fcs.h"

void rref(const int nrows,
          const int ncols,
          double **mat,
          int &nrank,
          const double tolerance = 1.0e-12);

void rref(std::vector<std::vector<double>> &mat,
          const double tolerance = 1.0e-12);

void rref_sparse(const size_t ncols,
                 ConstraintSparseForm &sp_constraint,
                 const double tolerance = 1.0e-12);
