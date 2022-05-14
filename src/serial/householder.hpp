#ifndef HOUSEHOLDER_HPP
#define HOUSEHOLDER_HPP

#include "../utils/data_structures.hpp"
#include "../helpers/helpers.hpp"
#include <cmath>

void matmul_serial(cse402project::matrix* Aptr, cse402project::matrix* Bptr, cse402project::matrix* resptr);

void householder_serial(cse402project::matrix* Aptr, cse402project::matrix* Hptr);

#endif