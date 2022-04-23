#include "utils/data_structures.hpp"
#include "utils/random_generators.hpp"

#include "serial/householder.hpp"

#include "openmp/householder_omp.hpp"

#include <omp.h>

#include <iomanip>


using cse402project::matrix;
using cse402project::vector;

int main(int argc, char **argv){

    //omp_set_num_threads(6);
    int n = 512;
    omp_set_num_threads(6);   
    matrix m1(n,n);

    double s1,t1,s2,t2;

    //Test data
    MATRIX_ELEMENT(m1,0,0) = 4;
    MATRIX_ELEMENT(m1,1,1) = 2;
    MATRIX_ELEMENT(m1,2,2) = 3;
    MATRIX_ELEMENT(m1,3,3) = -1;
    MATRIX_ELEMENT(m1,1,0) = 1;
    MATRIX_ELEMENT(m1,2,0) = -2;
    MATRIX_ELEMENT(m1,3,0) = 2;
    MATRIX_ELEMENT(m1,2,1) = 0;
    MATRIX_ELEMENT(m1,3,1) = 1;
    MATRIX_ELEMENT(m1,3,2) = -2;

    for(int i = 0; i<n; ++i){
        for(int j=i+1; j<n; ++j){
            MATRIX_ELEMENT(m1,i,j) = MATRIX_ELEMENT(m1,j,i);
        }
    }

    rand_symm_matrix_generator(&m1);
    //m1.print_matrix();
    matrix H1(n,n);
    matrix H2(n,n);
    //rand_symm_matrix_generator(&m1);

    s1 = omp_get_wtime();
    householder_serial(&m1,&H1);
    t1 = omp_get_wtime() - s1;

    s2 = omp_get_wtime();
    householder_omp(&m1,&H2);
    t2 = omp_get_wtime() - s2;

    //m1.print_matrix();
    //H1.print_matrix();
    //H2.print_matrix();

    std::cout<<t1<<std::endl;
    std::cout<<t2<<std::endl;
    return 0;
}