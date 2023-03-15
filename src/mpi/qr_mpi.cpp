#include "qr_mpi.hpp"
#include <iostream>
#include <float.h>

#define MAX_ITER 25

using cse402project::matrix;
using cse402project::vector;

double vec_l2norm_mpi(vector* vec){

    int my_rank, num_pes;


    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_pes);

    //Divide up the array almost equally. Scatterv will be used.
    int vecSize = vec->size;
    int extra_elems = vecSize%num_pes;
    int displ[num_pes];
    int send_counts[num_pes];
    int temp = 0;

    if(num_pes > vecSize){

        std::cout<<"Number of processes cannot be greater than vector size!"<<std::endl;
        exit(1);

    }

    for(int i = 0; i < num_pes; ++i){
        send_counts[i] = vecSize/num_pes;
        if(extra_elems > 0){
            send_counts[i]++;
            extra_elems--;
        }

        displ[i] = temp;
        temp+=send_counts[i];
    }

    int recv_count = send_counts[my_rank];

    double recvbuf[recv_count];

    MPI_Scatterv(vec->data,send_counts,displ,MPI_DOUBLE,recvbuf,recv_count,MPI_DOUBLE,0, MPI_COMM_WORLD);

    double l2norm, sum_local=0, sum_global=0;

    for(int i=0; i<recv_count; ++i){

        sum_local = sum_local + recvbuf[i]*recvbuf[i];
    }

    MPI_Reduce(&sum_local, &sum_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    l2norm = sqrt(sum_global);

    return l2norm;

}

void matmul_mpi(matrix* A, matrix* B, matrix* res){

    int my_rank, num_pes;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_pes);

    int matrixSize = A->rows;

    int sum = 0;
    int extra_rows = matrixSize%num_pes;
    int displ[num_pes];
    int send_counts[num_pes];
    int rows_proc[num_pes];

    if(num_pes > matrixSize){
        std::cout<<"Number of processes cannot be greater than number of rows!"<<std::endl;
        exit(1);
    }

    for(int i = 0; i < num_pes; ++i){
        send_counts[i] = matrixSize*(matrixSize/num_pes);
        rows_proc[i] = matrixSize/num_pes;
        if(extra_rows > 0){
            send_counts[i]+=matrixSize;
            rows_proc[i]++;
            extra_rows--;
        }

        displ[i] = sum;
        sum+=send_counts[i];
    }

    int recv_count = send_counts[my_rank];
    int num_rows = rows_proc[my_rank];

    std::vector<double> recvbuf_A(num_rows*matrixSize);
    std::vector<double> sendbuf_C;
    sendbuf_C.reserve(num_rows*matrixSize);

    //double recvbuf_A[num_rows*matrixSize];
    //double sendbuf_C[num_rows*matrixSize] = {0};

    //std::cout<<"Rank = "<<my_rank<<" Rows = "<<num_rows<<" sendcounts = "<<send_counts[my_rank]<<" disp= "<<displ[my_rank]<<std::endl;

    MPI_Scatterv(A->data, send_counts, displ, MPI_DOUBLE, &recvbuf_A[0], recv_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(B->data, matrixSize*matrixSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    double Aik;
    for(int i = 0; i < num_rows; ++i){
        for(int k = 0; k < matrixSize; ++k){

            Aik = recvbuf_A[k + i*matrixSize];
            for(int j = 0; j < matrixSize; ++j){
                sendbuf_C[j + matrixSize*i] += Aik*MATPTR_ELEMENT(B,k,j);
            }
        }
    }

    MPI_Gatherv(&sendbuf_C[0], num_rows*matrixSize, MPI_DOUBLE, res->data,send_counts, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
}

void householder_mpi(matrix* Aptr, matrix* Hptr){

    int my_rank, num_pes;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_pes);

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
        
        alpha = (double) -sgn(VECTOR_ELEMENT(colVec,k+1))*vec_l2norm_mpi(&colVec);

        if(my_rank == 0){
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
        }
        

        matmul_mpi(&P,&Ak,&temp);
        matmul_mpi(&temp, &P, &Ak);
    }

    *Hptr = Ak;

}

void qr_factorization_tridiagonal_mpi(matrix* A, matrix* Q, matrix* R){

    int my_rank, num_pes;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_pes);

    double x,y,r;
    double c,s;

    matrix Gi(A->rows,A->cols);
    matrix Gi_T(A->rows, A->cols);

    matrix temp(A->rows,A->cols);
    matrix tempQ(A->rows,A->cols);

    identity_matrix_serial(Q);

    matrix Ai(*A);

    for(int i=0; i<A->cols-1; ++i){

        if(my_rank == 0){
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
        }

        matmul_mpi(&Gi, &Ai, &temp);

        matmul_mpi(Q,&Gi_T,&tempQ);

        Ai = temp;
        *Q = tempQ;

    }

    *R = Ai;

}

void qr_iteration_mpi(matrix* A, vector* eigen_vals){

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

        qr_factorization_tridiagonal_mpi(&Ai, &Q, &R);
        matmul_mpi(&R,&Q,&Ai);

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


