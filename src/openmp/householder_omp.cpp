#include "householder_omp.hpp"
#include "../helpers/helpers.hpp"

using cse402project::matrix;
using cse402project::vector;

double vec_l2norm_omp(vector* vecptr){

    double l2norm = 0.0;
    double sum_square = 0.0;

    #pragma omp parallel for reduction(+:sum_square)
    for(int i=0; i<vecptr->size; ++i){

        sum_square = sum_square + VECPTR_ELEMENT(vecptr,i)*VECPTR_ELEMENT(vecptr,i);

    }

    l2norm = sqrt(sum_square);
    return l2norm;
}

void identity_matrix_omp(matrix* matptr){

    if(matptr->rows != matptr->cols){
        std::cout<<"Need a square matrix!"<<std::endl;
        exit(1);
    }

    int i;
    #pragma omp parallel for 
    for(i = 0; i<matptr->rows; ++i){
        MATPTR_ELEMENT(matptr,i,i) = 1.0;
    }

}

void matmul_tiled_omp(matrix* Aptr, matrix* Bptr, matrix* resptr, int tilesize){
    int i,j,k;
    int ii,jj,kk;

    //A temporary variable.
    double Aik;

    resptr->clear_matrix();

    #pragma omp parallel for private(kk,jj,i,j,k, Aik)
    for(ii=0; ii<Aptr->rows; ii+=tilesize){
        for(kk = 0; kk < Aptr->cols; kk+=tilesize){
            for(jj = 0; jj < Bptr->cols; jj+=tilesize){

                for(i = ii; i<ii+tilesize; ++i){
                    for(k = kk; k < kk + tilesize; ++k){

                        Aik = MATPTR_ELEMENT(Aptr,i,k);
                        for(j = jj; j < jj+tilesize; ++j){
                            MATPTR_ELEMENT(resptr,i,j) += Aik*MATPTR_ELEMENT(Bptr,k,j);
                        }

                    }
                }
            }
        }
    }
}