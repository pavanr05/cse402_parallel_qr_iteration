#include "mpi_manager.hpp"

void setup_env(int argc, char** argv, int &my_rank, int& num_pes){

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &num_pes);
  
}


