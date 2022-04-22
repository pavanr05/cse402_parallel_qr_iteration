#include "utils/data_structures.hpp"
#include "utils/random_generators.hpp"

#include "serial/householder.hpp"

#include "openmp/householder_omp.hpp"

#include <omp.h>

#include <iomanip>


using cse402project::matrix;
using cse402project::vector;

int main(int argc, char **argv)
{   
    matrix m1(3,3);
    rand_symm_matrix_generator(&m1);
    m1.print_matrix();
    return 0;
}