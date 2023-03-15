#ifndef QR_MPI_HPP
#define QR_MPI_HPP

#include "../utils/data_structures.hpp"
#include "../utils/utils.hpp"
#include "../helpers/helpers.hpp"
#include <vector>
#include <mpi.h>

//1D tiled matrix multiplication.
void matmul_mpi(cse402project::matrix* A, cse402project::matrix* B, cse402project::matrix* res);

void householder_mpi(cse402project::matrix* A, cse402project::matrix* H);

double vec_l2norm_mpi(cse402project::vector* vec);

void qr_factorization_tridiagonal_mpi(cse402project::matrix* A, cse402project::matrix* Q, cse402project::matrix* R);

void qr_iteration_mpi(cse402project::matrix* A, cse402project::vector* eigen_vals);
#endif