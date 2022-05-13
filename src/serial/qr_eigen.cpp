#include "qr_eigen.hpp"
#include <algorithm>
#include <cmath>
#include <float.h>

#define rel_tol 1e-1
#define MAX_ITER 25

using cse402project::matrix;
using cse402project::vector;

void qr_factorization_tridiagonal(matrix* A, matrix* Q, matrix* R){

    double x,y,r;
    double c,s;

    matrix Gi(A->rows,A->cols);
    matrix Gi_T(A->rows, A->cols);

    matrix temp(A->rows,A->cols);
    matrix tempQ(A->rows,A->cols);

    identity_matrix_serial(Q);

    matrix Ai(*A);

    for(int i=0; i<A->cols-1; ++i){

        identity_matrix_serial(&Gi);
        identity_matrix_serial(&Gi_T);

        x = MATRIX_ELEMENT(Ai,i,i);
        y = MATRIX_ELEMENT(Ai,i+1,i);
        r = sqrt(x*x + y*y);

        c = x/r;
        s = -y/r;

        MATRIX_ELEMENT(Gi,i,i) = c;
        MATRIX_ELEMENT(Gi,i,i+1) = -s;
        MATRIX_ELEMENT(Gi,i+1,i) = s;
        MATRIX_ELEMENT(Gi,i+1,i+1) = c;

        MATRIX_ELEMENT(Gi_T,i,i) = c;
        MATRIX_ELEMENT(Gi_T,i,i+1) = s;
        MATRIX_ELEMENT(Gi_T,i+1,i) = -s;
        MATRIX_ELEMENT(Gi_T,i+1,i+1) = c;

        matmul_serial(&Gi, &Ai, &temp);

        matmul_serial(Q,&Gi_T,&tempQ);

        Ai = temp;
        *Q = tempQ;

    }

    *R = Ai;
}


void qr_iteration(matrix* A, vector* eigen_vals){

    matrix Q(A->rows,A->cols);
    matrix R(A->rows,A->cols);
    matrix Ai(*A);

    int numIters = 0;
    double rel_err = DBL_MAX;
    double err_vec_norm;
    double v_prev_norm;

    vector v_prev(A->rows);
    std::fill_n(v_prev.data, v_prev.size, 1);
    vector err_vec(A->rows);


    while(numIters<MAX_ITER){

        qr_factorization_tridiagonal(&Ai, &Q, &R);
        matmul_serial(&R,&Q,&Ai);

        /*
        for(int i = 0; i<Ai.rows; ++i){
            VECTOR_ELEMENT(err_vec,i) = MATRIX_ELEMENT(Ai,i,i) - VECTOR_ELEMENT(v_prev,i);
        }

        err_vec_norm = vec_l2norm_serial(&err_vec);
        v_prev_norm = vec_l2norm_serial(&v_prev);
        
        rel_err = err_vec_norm/v_prev_norm;
        */
       
        for(int i = 0; i<A->rows; ++i){
            VECTOR_ELEMENT(v_prev,i) = MATRIX_ELEMENT(Ai,i,i);
        }
        

        numIters++;


    }

    *eigen_vals = v_prev;

}

