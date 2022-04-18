#ifndef DATA_STRUCTURES_HPP
#define DATA_STRUCTURES_HPP

#include <iostream>

namespace cse402project{

    class vector{

        public:
            
            double *data;
            int size;

            vector();

            //Construct a vector of size_ initialized with 0.
            vector(int size_);

            //Copy constructor
            vector(const vector& other);

            //Assignment operator
            vector& operator=(const vector& other);

            //Destructor
            ~vector();

            void print_vector();

            void clear_vector();

            private:

                //Helper to delete vector.
                void delete_vector();
    };

    class matrix{

        public:

            double **data;
            int rows;
            int cols;

            matrix();

            //Construct a matrix with rows_ and cols_ initialized to 0.
            matrix(int rows_, int cols_);

            //Copy constructor
            matrix(const matrix& other);

            //Assignment operator
            matrix& operator=(const matrix& other);

            //Destructor
            ~matrix();

            void print_matrix();

            //Initialize all elements to 0
            void clear_matrix();

            private:

                //Helper to create matrix.
                void create_matrix();

                //Helper to delete matrix.
                void delete_matrix();
    };
}

#define VECPTR_ELEMENT(vecptr,i) vecptr->data[i]
#define VECTOR_ELEMENT(vec, i) vec.data[i]

#define MATPTR_ELEMENT(matptr, i, j) matptr->data[i][j]
#define MATRIX_ELEMENT(mat, i, j) mat.data[i][j]

#endif