#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>

namespace cse402project{

    struct vector{
        double *data;
        int size;
    };

    struct matrix{
        double *data;
        int rows;
        int cols;
    };
}

#define ARRAY_ELEMENT(dptr, N, M, i, j) dptr[(i)*(M) + (j)]
#define MATPTR_ELEMENT(matptr, i, j) ARRAY_ELEMENT(matptr->data, matptr->rows, matptr->cols, i, j)
#define MATRIX_ELEMENT(matr, i, j) ARRAY_ELEMENT(matr.data, matr.rows, matr.cols, i, j)
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

bool matrix_create(cse402project::matrix *target, int rows, int cols);
bool vector_create(cse402project::vector *target, int size);

bool matrix_destroy(cse402project::matrix *target);
bool vector_destroy(cse402project::vector *target);

bool copy_matrix(cse402project::matrix *target, cse402project::matrix *source);
bool copy_vector(cse402project::vector *target, cse402project::vector *source);


#endif