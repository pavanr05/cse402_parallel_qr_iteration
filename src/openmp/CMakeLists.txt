cmake_minimum_required (VERSION 3.5)

find_package(OpenMP)

add_library(helpersomp STATIC ../helpers/helpers.cpp ../helpers/helpers.hpp)


add_library(utilsomp STATIC ../utils/data_structures.cpp ../utils/data_structures.hpp ../utils/utils.cpp ../utils/utils.hpp)

add_executable(eigen_omp main_omp.cpp householder_omp.cpp householder_omp.hpp)

target_link_libraries(utilsomp helpersomp)


if(OpenMP_CXX_FOUND)
    target_link_libraries(eigen_omp PUBLIC utilsomp OpenMP::OpenMP_CXX)
endif()

include_directories(../helpers ../utils)

target_compile_options(eigen_omp PUBLIC -qopenmp)