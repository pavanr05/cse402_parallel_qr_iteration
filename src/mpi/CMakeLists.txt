cmake_minimum_required (VERSION 3.5)

string(FIND "${CMAKE_CXX_COMPILER}" "mpiicpc" MPIICPC)
string(FIND "${CMAKE_CXX_COMPILER}" "mpicxx" MPICXX)
string(FIND "${CMAKE_CXX_COMPILER}" "mpic++" MPICPP)

math(EXPR MPIICPC_INT ${MPIICPC})
math(EXPR MPICXX_INT ${MPICXX}) 
math(EXPR MPICPP_INT ${MPICPP})


message(STATUS "${CMAKE_CXX_COMPILER}")
find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})
if((MPIICPC_INT GREATER_EQUAL 0) OR (MPICXX_INT GREATER_EQUAL 0) OR (MPICPP_INT GREATER_EQUAL 0))

    target_link_libraries(utils ${MPI_CXX_LIBRARIES})

    if(MPI_COMPILE_FLAGS)
    set_target_properties(utils PROPERTIES
    COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    endif()


    add_executable(eigen_mpi main_mpi.cpp mpi_manager.hpp mpi_manager.cpp qr_mpi.hpp qr_mpi.cpp)

    target_link_libraries(eigen_mpi utils)
endif()
