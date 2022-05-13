#ifndef QR_MPI_HPP
#define QR_MPI_HPP

#include "../utils/data_structures.hpp"
#include "../utils/utils.hpp"
#include <mpi.h>


//1D tiled matrix multiplication.
void matmul_mpi(cse402project::matrix* A, cse402project::matrix* B, cse402project::matrix* res);

#endif