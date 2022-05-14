#include "householder.hpp"

using cse402project::matrix;
using cse402project::vector;


void matmul_serial(matrix* Aptr, matrix* Bptr, matrix* resptr){
    
    //A temporary variable.
    double Aik;

    resptr->clear_matrix();

    for(int k = 0; k<Aptr->cols; ++k){
        for(int i = 0; i < Aptr->rows; ++i){
            Aik = MATPTR_ELEMENT(Aptr,i,k);
            for(int j=0; j < Bptr->cols; ++j){
                MATPTR_ELEMENT(resptr,i,j) += Aik*MATPTR_ELEMENT(Bptr,k,j);
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

        matmul_serial(&P,&Ak,&temp);
        matmul_serial(&temp, &P, &Ak);
    }

    *Hptr = Ak;

}

