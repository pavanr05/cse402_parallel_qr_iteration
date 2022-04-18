#include "data_structures.hpp"

int main(){

    cse402project::matrix m1(10,10);
    MATRIX_ELEMENT(m1,2,2) = 5;
    cse402project::matrix m2;
    m2 = m1;
    m1.clear_matrix();
    m1.print_matrix();
    m2.print_matrix();
}