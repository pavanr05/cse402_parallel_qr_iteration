#ifndef HOUSEHOLDER_OMP_HPP
#define HOUSEHOLDER_OMP_HPP

#include "../utils/data_structures.hpp"
#include <cmath>
#include <omp.h>



//Calculate the L2 norm of the vector.
double vec_l2norm_omp(cse402project::vector* vecptr);

#endif