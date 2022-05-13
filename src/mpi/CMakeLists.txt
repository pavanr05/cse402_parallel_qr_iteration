cmake_minimum_required (VERSION 3.5)

find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

target_link_libraries(utils ${MPI_CXX_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(utils PROPERTIES
  COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()


add_executable(eigen_mpi main_mpi.cpp mpi_manager.hpp mpi_manager.cpp qr_mpi.hpp qr_mpi.cpp)

target_link_libraries(eigen_mpi utils)
