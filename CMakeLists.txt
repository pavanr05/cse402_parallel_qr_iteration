cmake_minimum_required(VERSION 3.5)
project(parallel_eigen LANGUAGE CXX)

add_subdirectory (src/mpi)
add_subdirectory (src/openmp)
add_subdirectory (src/serial)
add_subdirectory (tests)