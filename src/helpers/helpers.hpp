#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "../utils/data_structures.hpp"
#include "../utils/utils.hpp"
#include <cmath>


int sgn(double a);

//Calculate the L2 norm of the vector.
double vec_l2norm_serial(cse402project::vector* vecptr);

void identity_matrix_serial(cse402project::matrix* matptr);

#endif