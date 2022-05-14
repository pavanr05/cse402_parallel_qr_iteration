#include "helpers.hpp"
#include <iostream>

using cse402project::matrix;
using cse402project::vector;

int sgn(double a){
    if(a>=0) return 1;
    else return -1;
}

double vec_l2norm_serial(vector* vecptr){

    double l2norm = 0.0;
    double sum_square = 0.0;

    for(int i=0; i<vecptr->size; ++i){

        sum_square = sum_square + VECPTR_ELEMENT(vecptr,i)*VECPTR_ELEMENT(vecptr,i);

    }

    l2norm = sqrt(sum_square);
    return l2norm;
}

void identity_matrix_serial(matrix* matptr){

    if(matptr->rows != matptr->cols){
        std::cout<<"Need a square matrix!"<<std::endl;
        exit(1);
    }

    matptr->clear_matrix();
    
    for(int i = 0; i<matptr->rows; ++i){
        MATPTR_ELEMENT(matptr,i,i) = 1.0;
    }
}