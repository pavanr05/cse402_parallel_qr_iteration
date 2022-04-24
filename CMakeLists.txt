cmake_minimum_required(VERSION 3.5)
set(CMAKE_C_COMPILER mpiicc)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_CXX_STANDARD 11)

project(eigen_parallel)

add_subdirectory(src/serial)
add_subdirectory(src/openmp)

