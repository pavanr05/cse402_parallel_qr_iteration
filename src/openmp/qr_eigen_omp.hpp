#ifndef QR_EIGEN_OMP
#define QR_EIGEN_OMP

#include "../helpers/helpers.hpp"
#include "../utils/data_structures.hpp"
#include "../utils/utils.hpp"
#include "householder_omp.hpp"
#include <omp.h>

//Note that this function expects a tridiagonal matrix A.
void qr_factorization_tridiagonal_omp(cse402project::matrix* A, cse402project::matrix* Q, cse402project::matrix* R);

//QR iteration to find the eigen values.
void qr_iteration_omp(cse402project::matrix* A, cse402project::vector* eigen_vals);
#endif