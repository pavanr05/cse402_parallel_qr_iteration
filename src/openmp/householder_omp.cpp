#include "householder_omp.hpp"

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

    matptr->clear_matrix();
    int i;
    #pragma omp parallel for 
    for(i = 0; i<matptr->rows; ++i){
        MATPTR_ELEMENT(matptr,i,i) = 1.0;
    }

}

void matmul_omp(matrix* A, matrix* B, matrix* res){
    double Aik;

    int i,j,k;

    res->clear_matrix();

    #pragma omp parallel for private(j,k,Aik)
    for(i = 0; i < A->rows; ++i){
        for(k = 0; k < A->cols; ++k){

            Aik = MATPTR_ELEMENT(A,i,k);
            for(j = 0; j < B->cols; ++j){
                MATPTR_ELEMENT(res,i,j) += Aik*MATPTR_ELEMENT(B,k,j);
            }
        }
    }
}

void householder_omp(matrix* Aptr, matrix* Hptr){

    double alpha, r;
    vector v(Aptr->rows);
    vector colVec(Aptr->rows);
    
    matrix P(Aptr->rows, Aptr->cols);

    //Identity matrix
    matrix I(Aptr->rows, Aptr->cols);
    identity_matrix_omp(&I);

    matrix temp(Aptr->rows, Aptr->cols);

    matrix Ak(*Aptr);

    //This loop cannot be parallelized because of the dependency of the Ak matrix from previous iteration.
    for(int k = 0; k < Aptr->cols-2; ++k){

        identity_matrix_omp(&P);

        std::fill_n(colVec.data, k+1, 0);

        #pragma omp parallel for
        for(int j = k+1; j < colVec.size; ++j){
            VECTOR_ELEMENT(colVec,j) = MATRIX_ELEMENT(Ak,j,k);
        }
        
        alpha = (double) -sgn(VECTOR_ELEMENT(colVec,k+1))*vec_l2norm_omp(&colVec);

        r = sqrt(0.5*alpha*alpha - 0.5*VECTOR_ELEMENT(colVec,k+1)*alpha);

        std::fill_n(v.data,k+1,0);

        VECTOR_ELEMENT(v,k+1) = 0.5*(VECTOR_ELEMENT(colVec,k+1) - alpha)/r;

        #pragma omp parallel for
        for(int j = k+2; j < Aptr->rows; ++j){
            VECTOR_ELEMENT(v,j) = 0.5*VECTOR_ELEMENT(colVec,j)/r;
        }

        #pragma omp parallel for
        for(int i=k+1; i < Aptr->rows; ++i){
            for(int j = k+1; j < Aptr->rows; ++j){
                MATRIX_ELEMENT(P,i,j) = MATRIX_ELEMENT(I,i,j) - 2*VECTOR_ELEMENT(v,i)*VECTOR_ELEMENT(v,j);
            }
        }

        matmul_omp(&P,&Ak,&temp);
        matmul_omp(&temp, &P, &Ak);
    }

    *Hptr = Ak;

}