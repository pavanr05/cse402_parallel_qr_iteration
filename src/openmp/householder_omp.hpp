#ifndef HOUSEHOLDER_OMP_HPP
#define HOUSEHOLDER_OMP_HPP

#include "../utils/data_structures.hpp"
#include "../helpers/helpers.hpp"
#include <cmath>
#include <omp.h>



//Calculate the L2 norm of the vector.
double vec_l2norm_omp(cse402project::vector* vecptr);

void identity_matrix_omp(cse402project::matrix*matptr);

void matmul_omp(cse402project::matrix* A, cse402project::matrix* B, cse402project::matrix* res);

void householder_omp(cse402project::matrix* Aptr, cse402project::matrix* Hptr);

#endif