#ifndef QR_EIGEN_HPP
#define QR_EIGEN_HPP

#include "../utils/data_structures.hpp"
#include "../helpers/helpers.hpp"
#include "householder.hpp"

//Note that this function expects a tridiagonal matrix A.
void qr_factorization_tridiagonal(cse402project::matrix* A, cse402project::matrix* Q, cse402project::matrix* R);

//Untiled matrix multiplication for small matrices.
void matmul(cse402project::matrix* Aptr, cse402project::matrix* Bptr, cse402project::matrix* resptr);

//QR iteration to find the eigen values.
void qr_iteration(cse402project::matrix* A, cse402project::vector* eigen_vals);
#endif