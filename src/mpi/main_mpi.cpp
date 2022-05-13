#include "mpi_manager.hpp"
#include "../utils/data_structures.hpp"
#include "../utils/utils.hpp"
#include "qr_mpi.hpp"

#include <stdlib.h>
#include <mpi.h>
#include <thread>
#include <iostream>

using cse402project::matrix;
using cse402project::vector;

int main(int argc, char** argv){

    if(argc<2){
        std::cout<<"Enter the matrix size"<<std::endl;
        exit(1);
    }

    int matrixSize = atoi(argv[1]);

    int my_rank, num_pes;

    setup_env(argc, argv, my_rank, num_pes);

    matrix m1(matrixSize, matrixSize);
    matrix m2(matrixSize, matrixSize);
    matrix m3(matrixSize, matrixSize);

    rand_symm_matrix_generator(&m1);
    rand_matrix_generator(&m2);
    matmul_mpi(&m1,&m2, &m3);
    

    MPI_Finalize();
    return 0;
}