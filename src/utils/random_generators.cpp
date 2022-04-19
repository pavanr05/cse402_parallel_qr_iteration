#include "random_generators.hpp"

void rand_matrix_generator(matrix *matptr){
    int randdivide = rand()%10000;
    
    for(int i=0; i<matptr->rows; ++i){
        for(int j=0; j<matptr->cols; ++j){
            MATPTR_ELEMENT(matptr,i,j) = (double) rand() / (double) randdivide;
        }
    }
}

void rand_vector_generator(vector* vecptr){
    int randdivide = rand()%10000;

    for(int i=0; i<vecptr->size; ++i){
        VECPTR_ELEMENT(vecptr,i) = (double) rand() / (double) randdivide;
    }

}

void rand_symm_matrix_generator(matrix* matptr){

    if(matptr->rows != matptr->cols){
        std::cout<<"Cannot generate a non-square symmteric matrix!"<<std::endl;
        exit(1);
    }

    int randdivide = rand()%10000;

    for(int i=0; i<matptr->rows; ++i){
        for(int j=0; j<=i; ++j){
            MATPTR_ELEMENT(matptr,i,j) = (double) rand() / (double) randdivide;
        }
    }

    for(int i=0; i<matptr->rows; ++i){
        for(int j=i+1; j<matptr->cols; ++j){
            MATPTR_ELEMENT(matptr,i,j) = MATPTR_ELEMENT(matptr,j,i);
        }
    }

}