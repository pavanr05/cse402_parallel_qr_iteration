#include "data_structures.hpp"
#include <unistd.h>

using cse402project::matrix;
using cse402project::vector;

//Taken from CS 420 MP1.
bool matrix_create(matrix *target, int rows, int cols){
    target->rows = rows;
    target->cols = cols;

    if(rows==0 || cols==0){
        target->data = NULL;
    }

    int memalign_result = posix_memalign(
            (void**)(&(target->data)),
            sysconf(_SC_LEVEL1_DCACHE_LINESIZE),
            sizeof(double)*rows*cols
    );
    if(NULL == target->data){
        target->rows = 0;
        target->cols = 0;
        return false;
    }

    target->rows = rows;
    target->cols = cols;
    return true;
}

bool matrix_destroy(matrix *target){
    free(target->data);
    target->data = NULL;
    target->rows = 0;
    target->cols = 0;
    return true;
}

bool copy_matrix(matrix *target, matrix *source){
    if(NULL == target->data || NULL == source->data){
        return false;
    }

    if(target->rows*target->cols != source->rows*source->cols){
        matrix_destroy(target);
        bool recreate = matrix_create(target,source->rows, source->cols);
        if(!recreate)
            return false;
    }

    target->rows = source->rows;
    target->cols = source->cols;

    for(int i = 0;i<source->rows;i++)
        for(int j = 0;j<source->cols;j++)
            MATPTR_ELEMENT(target, i, j) = MATPTR_ELEMENT(source, i, j);

    return true;
}

