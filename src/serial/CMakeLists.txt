cmake_minimum_required (VERSION 3.5)

add_executable(eigen_serial main_serial.cpp householder.cpp householder.hpp)

target_link_libraries(eigen_serial utils)

