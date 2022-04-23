#include "householder.hpp"
#include "../utils/data_structures.hpp"
#include "../utils/utils.hpp"
#include <stdlib.h>

using cse402project::matrix;
using cse402project::vector;

int main(int argc, char **argv){

    if(argc<2){
        std::cout<<"Enter the matrix size!"<<std::endl;
        exit(1);
    }

    int matrixSize = atoi(argv[1]);
    
    matrix m1(matrixSize,matrixSize);
    matrix H1(matrixSize,matrixSize);

    rand_symm_matrix_generator(&m1);

    double start,exec_time;

    start = get_wall_time();
    householder_serial(&m1, &H1);
    exec_time = get_wall_time() - start;

    std::cout<<"The serial version took "<<exec_time<<" s."<<std::endl;


    return 0;
}