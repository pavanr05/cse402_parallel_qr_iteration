#include "householder_omp.hpp"
#include "../utils/data_structures.hpp"
#include "../utils/utils.hpp"
#include "qr_eigen_omp.hpp"
#include <stdlib.h>
#include <omp.h>
#include <thread>

using cse402project::matrix;
using cse402project::vector;

int main(int argc, char **argv){

    if(argc<2){
        std::cout<<"Enter the matrix size"<<std::endl;
        exit(1);
    }

    const auto num_cores = std::thread::hardware_concurrency();

    int matrixSize = atoi(argv[1]);
    int numThreads;

    if(argc==3){
        numThreads = atoi(argv[2]);
    }else if(argc==2){
        numThreads = num_cores/2;
    }else{
        std::cout<<"Invalid arguments!"<<std::endl;
        exit(1);
    }

    omp_set_num_threads(numThreads);
    

    matrix m1(matrixSize,matrixSize);
    matrix H1(matrixSize,matrixSize);

    matrix Q(matrixSize, matrixSize);
    matrix R(matrixSize, matrixSize);

    vector eigen_vals(matrixSize);

    rand_symm_matrix_generator(&m1);
    //rand_matrix_generator(&Q);

    /*
    matmul_omp(&m1, &Q, &R);

    m1.print_matrix();
    Q.print_matrix();
    R.print_matrix();
    */

    
    //householder_omp(&m1, &H1);
    //H1.print_matrix();
    /*
    MATRIX_ELEMENT(H1,0,0) = 3;
    MATRIX_ELEMENT(H1,1,1) = 3;
    MATRIX_ELEMENT(H1,2,2) = 3;

    MATRIX_ELEMENT(H1,0,1) = 1;
    MATRIX_ELEMENT(H1,1,0) = 1;

    MATRIX_ELEMENT(H1,1,2) = 1;
    MATRIX_ELEMENT(H1,2,1) = 1;

    H1.print_matrix();

    qr_iteration_omp(&H1, &eigen_vals);

    eigen_vals.print_vector();
    */


    double start,exec_time;

    start = get_wall_time();
    householder_omp(&m1, &H1);
    qr_iteration_omp(&H1, &eigen_vals);
    exec_time = get_wall_time() - start;

    //eigen_vals.print_vector();

    std::cout<<"The OpenMP version with "<<numThreads<<" threads took "<<exec_time<<"s for a matrix of size "<<matrixSize<<"."<<std::endl;


    return 0;
}