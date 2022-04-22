#include "data_structures.hpp"
#include <algorithm>
#include <iomanip>

using cse402project::matrix;
using cse402project::vector;

//Vector class
void vector::create_vector(){
    data = new double[size];
}
void vector::clear_vector(){
    std::fill_n(data,size,0);
}

void vector::delete_vector(){

    delete [] data;
    size = 0;
}

vector::vector(){
    size = 0;
    data = NULL;
}

vector::vector(int size_){

    size = size_;
    create_vector();

    std::fill_n(data, size_, 0);
}

vector::vector(const vector& other){
    size = other.size;

    create_vector();
    for(int i = 0; i<size; ++i){
        data[i] = other.data[i];
    }
}

vector& vector::operator=(const vector& other){

    if(this != &other){

        delete_vector();
        size = other.size;

        create_vector();
        for(int i = 0; i<size; ++i){
            data[i] = other.data[i];
        }
    }
    
    return *this;
}

vector::~vector(){
    delete_vector();
}

void vector::print_vector(){
    for(int i=0; i < size; ++i){
        std::cout<<std::setprecision(PRECISION_DIGITS)<<data[i]<<" ";
    }

    std::cout<<std::endl;
}

//Matrix class

void matrix::create_matrix(){
    
    data = new double*[rows];
    for(int i = 0; i< rows; ++i){
        data[i] = new double[cols];
    }

}

void matrix::delete_matrix(){
    for(int i=0; i<rows; ++i){
        delete [] data[i];
    }

    rows = 0;
    cols = 0;
}

void matrix::clear_matrix(){
    for(int i=0; i<rows; ++i){
        std::fill_n(data[i],cols,0);
    }
}

matrix::matrix(){
    rows = 0;
    cols = 0;
    data = NULL;
}

matrix::matrix(int rows_, int cols_){

    rows = rows_;
    cols = cols_;
    create_matrix();
    clear_matrix();
}

matrix::matrix(const matrix& other){

    rows = other.rows;
    cols = other.cols;
    create_matrix();

    for(int i=0; i<rows; ++i){
        for(int j=0; j<cols; ++j){
            data[i][j] = other.data[i][j];
        }
    }
}

matrix& matrix::operator=(const matrix& other){
    if(this != &other){

        delete_matrix();

        rows = other.rows;
        cols = other.cols;

        create_matrix();

        for(int i=0; i<rows; ++i){
            for(int j=0; j<cols; ++j){
                data[i][j] = other.data[i][j];
            }
        }
    }

    return *this;
}

matrix::~matrix(){
    delete_matrix();
}

void matrix::print_matrix(){

    for(int i=0; i<rows; ++i){
        for(int j=0; j<cols; ++j){
            std::cout<<std::setprecision(PRECISION_DIGITS)<<data[i][j]<<" ";
        }

        std::cout<<std::endl;
    }

    std::cout<<std::endl;
}