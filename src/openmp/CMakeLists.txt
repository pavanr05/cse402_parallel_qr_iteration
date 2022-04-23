cmake_minimum_required (VERSION 3.5)
set (CMAKE_C_COMPILER gcc)
set (CMAKE_C_STANDARD 11)
set (CMAKE_CXX_COMPILER g++)
set (CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall -O3)


add_library(helpersomp STATIC ../helpers/helpers.cpp ../helpers/helpers.hpp)


add_library(utilsomp STATIC ../utils/data_structures.cpp ../utils/data_structures.hpp ../utils/utils.cpp ../utils/utils.hpp)

add_executable(eigen_omp main_omp.cpp householder_omp.cpp householder_omp.hpp)

target_link_libraries(utilsomp helpersomp)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(eigen_omp PUBLIC utilsomp OpenMP::OpenMP_CXX)
endif()

include_directories(../helpers ../utils)

target_compile_options(eigen_omp PUBLIC -fopenmp)