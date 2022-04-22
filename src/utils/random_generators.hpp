#ifndef RANDOM_GENERATORS_HPP
#define RANDOM_GENERATORS_HPP

#include "data_structures.hpp"
#include <cstdlib>



//Random matrix generator.
void rand_matrix_generator(cse402project::matrix *matptr);

//Random vector generator.
void rand_vector_generator(cse402project::vector *vecptr);

//Random symmetric matrix.
void rand_symm_matrix_generator(cse402project::matrix *matptr);

#endif