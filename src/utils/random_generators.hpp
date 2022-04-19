#ifndef RANDOM_GENERATORS_HPP
#define RANDOM_GENERATORS_HPP

#include "data_structures.hpp"
#include <cstdlib>

using cse402project::matrix;
using cse402project::vector;

//Random matrix generator.
void rand_matrix_generator(matrix *matptr);

//Random vector generator.
void rand_vector_generator(vector *vecptr);

//Random symmetric matrix.
void rand_symm_matrix_generator(matrix *matptr);

#endif