#include "householder.hpp"
#include "qr_eigen.hpp"
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

    matrix Q(matrixSize, matrixSize);
    matrix R(matrixSize, matrixSize);

    vector eigen_vals(matrixSize);

    rand_symm_matrix_generator(&m1);

    /*
    MATRIX_ELEMENT(H1,0,0) = 3;
    MATRIX_ELEMENT(H1,1,1) = 3;
    MATRIX_ELEMENT(H1,2,2) = 3;

    MATRIX_ELEMENT(H1,0,1) = 1;
    MATRIX_ELEMENT(H1,1,0) = 1;

    MATRIX_ELEMENT(H1,1,2) = 1;
    MATRIX_ELEMENT(H1,2,1) = 1;
    */

    double start,exec_time;

    start = get_wall_time();
    householder_serial(&m1, &H1);
    qr_iteration(&H1, &eigen_vals);
    exec_time = get_wall_time() - start;

    //H1.print_matrix();
    //m1.print_matrix();
    //eigen_vals.print_vector();
    

    std::cout<<"The serial version took "<<exec_time<<"s for a matrix of size "<<matrixSize<<"."<<std::endl;

    


    return 0;
}