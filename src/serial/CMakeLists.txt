cmake_minimum_required (VERSION 3.5)

add_executable(eigen_serial main_serial.cpp householder.cpp householder.hpp qr_eigen.hpp qr_eigen.cpp)

target_link_libraries(eigen_serial utils)

