cmake_minimum_required(VERSION 3.5)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

project(eigen_parallel LANGUAGES CXX)

add_compile_options(-Wall -O3 -g)

add_library(helpers STATIC src/helpers/helpers.cpp src/helpers/helpers.hpp)
add_library(utils STATIC src/utils/data_structures.cpp src/utils/data_structures.hpp src/utils/utils.cpp src/utils/utils.hpp)

target_link_libraries(utils helpers)

include_directories(src/helpers src/utils)

add_subdirectory(src/serial)
add_subdirectory(src/openmp)

