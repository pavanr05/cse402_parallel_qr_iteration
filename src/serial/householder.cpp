#include "householder.hpp"
#define TILE_SIZE 8

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

void matmul_tiled(matrix* Aptr, matrix* Bptr, matrix* resptr, int tilesize){
    int i,j,k;
    int ii,jj,kk;

    //A temporary variable.
    double Aik;

    resptr->clear_matrix();

    for(kk=0; kk<Aptr->cols; kk+=tilesize){
        for(ii = 0; ii < Aptr->rows; ii+=tilesize){
            for(jj = 0; jj < Bptr->cols; jj+=tilesize){
                for(k = kk; k<kk+tilesize; ++k){
                    for(i = ii; i < ii + tilesize; ++i){

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

void householder_serial(matrix* Aptr, matrix* Hptr){

    double alpha, r;
    vector v(Aptr->rows);
    vector colVec(Aptr->rows);
    
    matrix P(Aptr->rows, Aptr->cols);

    //Identity matrix
    matrix I(Aptr->rows, Aptr->cols);
    identity_matrix_serial(&I);

    matrix temp(Aptr->rows, Aptr->cols);

    matrix Ak(*Aptr);

    for(int k = 0; k < Aptr->cols-2; ++k){

        identity_matrix_serial(&P);

        std::fill_n(colVec.data, k+1, 0);
        for(int j = k+1; j < colVec.size; ++j){
            VECTOR_ELEMENT(colVec,j) = MATRIX_ELEMENT(Ak,j,k);
        }
        
        alpha = (double) -sgn(VECTOR_ELEMENT(colVec,k+1))*vec_l2norm_serial(&colVec);

        r = sqrt(0.5*alpha*alpha - 0.5*VECTOR_ELEMENT(colVec,k+1)*alpha);

        std::fill_n(v.data,k+1,0);

        VECTOR_ELEMENT(v,k+1) = 0.5*(VECTOR_ELEMENT(colVec,k+1) - alpha)/r;

        for(int j = k+2; j < Aptr->rows; ++j){
            VECTOR_ELEMENT(v,j) = 0.5*VECTOR_ELEMENT(colVec,j)/r;
        }

        for(int i=k+1; i < Aptr->rows; ++i){
            for(int j = k+1; j < Aptr->rows; ++j){
                MATRIX_ELEMENT(P,i,j) = MATRIX_ELEMENT(I,i,j) - 2*VECTOR_ELEMENT(v,i)*VECTOR_ELEMENT(v,j);
            }
        }

        matmul_tiled(&P,&Ak,&temp,TILE_SIZE);
        matmul_tiled(&temp, &P, &Ak, TILE_SIZE);
    }

    *Hptr = Ak;

}

