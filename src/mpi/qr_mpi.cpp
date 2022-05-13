#include "qr_mpi.hpp"
#include <iostream>

using cse402project::matrix;
using cse402project::vector;




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
        std::cout<<"Number of processes cannt be greater than number of rows!"<<std::endl;
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

    double recvbuf_A[num_rows*matrixSize];
    double sendbuf_C[num_rows*matrixSize] = {0};

    std::cout<<"Rank = "<<my_rank<<" Rows = "<<num_rows<<" sendcounts = "<<send_counts[my_rank]<<" disp= "<<displ[my_rank]<<std::endl;

    MPI_Scatterv(A->data, send_counts, displ, MPI_DOUBLE, recvbuf_A, recv_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(B->data, matrixSize*matrixSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    std::cout<<std::endl;

    double Aik;
    for(int i = 0; i < num_rows; ++i){
        for(int k = 0; k < matrixSize; ++k){

            Aik = recvbuf_A[k + i*matrixSize];
            for(int j = 0; j < matrixSize; ++j){
                sendbuf_C[j + matrixSize*i] += Aik*MATPTR_ELEMENT(B,k,j);
            }
        }
    }

    MPI_Gatherv(sendbuf_C, num_rows*matrixSize, MPI_DOUBLE, res->data,send_counts, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
}