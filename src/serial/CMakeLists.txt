cmake_minimum_required (VERSION 3.5)

add_library(helpers STATIC ../helpers/helpers.cpp ../helpers/helpers.hpp)
add_library(utils STATIC ../utils/data_structures.cpp ../utils/data_structures.hpp ../utils/utils.cpp ../utils/utils.hpp)
add_executable(eigen_serial main_serial.cpp householder.cpp householder.hpp)

target_link_libraries(utils helpers)
target_link_libraries(eigen_serial utils)

include_directories(../helpers ../utils)