#include "householder_omp.hpp"
#include "../utils/data_structures.hpp"
#include "../utils/utils.hpp"
#include "qr_eigen_omp.hpp"
#include <stdlib.h>
#include <omp.h>

using cse402project::matrix;
using cse402project::vector;

int main(int argc, char **argv){

    if(argc<3){
        std::cout<<"Enter the matrix size and number of threads!"<<std::endl;
        exit(1);
    }

    int matrixSize = atoi(argv[1]);
    int numThreads = atoi(argv[2]);

    omp_set_num_threads(numThreads);
    
    matrix m1(matrixSize,matrixSize);
    matrix H1(matrixSize,matrixSize);

    matrix Q(matrixSize, matrixSize);
    matrix R(matrixSize, matrixSize);

    vector eigen_vals(matrixSize);

    rand_symm_matrix_generator(&m1);

    double start,exec_time;

    start = get_wall_time();
    householder_omp(&m1, &H1);
    qr_iteration_omp(&H1, &eigen_vals);
    exec_time = get_wall_time() - start;

    eigen_vals.print_vector();

    std::cout<<"The OpenMP version took "<<exec_time<<" s."<<std::endl;



    return 0;
}