#ifndef HOUSEHOLDER_HPP
#define HOUSEHOLDER_HPP

#include "../utils/data_structures.hpp"
#include "../helpers/helpers.hpp"
#include <cmath>



//Calculate the L2 norm of the vector.
double vec_l2norm_serial(cse402project::vector* vecptr);

void identity_matrix_serial(cse402project::matrix* matptr);

void matmul_tiled(cse402project::matrix* Aptr, cse402project::matrix* Bptr, cse402project::matrix* resptr, int tilesize);

void transpose_tiled(cse402project::matrix* Aptr, int tilesize);

void householder_serial(cse402project::matrix* Aptr, cse402project::matrix* Hptr);

#endif