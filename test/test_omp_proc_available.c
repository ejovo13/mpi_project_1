#include "geodesy.h"
#include "omp.h"
#include "mpi.h"

// Let's see how many threads are available on this machine

#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv){

    int this_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

    #pragma omp parallel
    {

        if (omp_get_thread_num() == 0) {
            printf("Node %d has available threads: %d\n", this_rank, omp_get_num_threads());
        }

    }

    MPI_Finalize();

    return 0;
}