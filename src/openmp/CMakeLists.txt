cmake_minimum_required (VERSION 3.5)

find_package(OpenMP)

add_executable(eigen_omp main_omp.cpp householder_omp.cpp householder_omp.hpp qr_eigen_omp.hpp qr_eigen_omp.cpp)

if(OpenMP_CXX_FOUND)
    target_link_libraries(eigen_omp PUBLIC utils OpenMP::OpenMP_CXX)
endif()

target_compile_options(eigen_omp PUBLIC -fopenmp)