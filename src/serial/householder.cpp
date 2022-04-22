#include "householder.hpp"

using cse402project::matrix;
using cse402project::vector;

double vec_l2norm_serial(vector* vecptr){

    double l2norm = 0.0;
    double sum_square = 0.0;

    for(int i=0; i<vecptr->size; ++i){

        sum_square = sum_square + VECPTR_ELEMENT(vecptr,i)*VECPTR_ELEMENT(vecptr,i);

    }

    l2norm = sqrt(sum_square);
    return l2norm;
}