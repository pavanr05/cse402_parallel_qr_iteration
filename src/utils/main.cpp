#include "data_structures.hpp"
#include "random_generators.hpp"
#include <algorithm>

using cse402project::matrix;
using cse402project::vector;

int main(){

    matrix m1(2,2);
    rand_symm_matrix_generator(&m1);
    m1.print_matrix();

    vector v1(10);
    rand_vector_generator(&v1);
    v1.print_vector();

    return 0;

    
}