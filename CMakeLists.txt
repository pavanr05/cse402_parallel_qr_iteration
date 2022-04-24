cmake_minimum_required(VERSION 3.5)

set(CMAKE_C_STARNDARD 11)
set(CMAKE_CXX_STANDARD 11)

project(parallel_eigen)

add_compile_options(-Wall -O3)

add_subdirectory(src/serial)
add_subdirectory(src/openmp)



