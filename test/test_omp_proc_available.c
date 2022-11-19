#include "geodesy.h"
#include "omp.h"
#include "mpi.h"

// Let's see how many threads are available on this machine

#include <stdio.h>
#include <mpi.h>

double my_func(double x) {
    return 3 * x * (x - 1);
}

int main(int argc, char** argv){

    int this_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

    // Create a matrix the size of the number of ranks (for gathering)
    Matrix_d *received_sums = Matrix_new_d(1, world_size);

    // Let's compute the sum of a function and then we'll sum that in a loop
    int N = 100000000;
    double sum = 0;
    double total_sum = 0;

    // Break up the sum based on this node's rank
    int workload = compute_workload(N, world_size, this_rank);
    Matrix_i *startend = compute_startend_array(N, world_size);
    int start = vecat_i(startend, this_rank) + 1;
    int end = vecat_i(startend, this_rank + 1);
    assert(workload == (end - start) + 1);
    if (this_rank == 0) {
        printf("Workloads asserted\n");
        Matrix_print_i(startend);
    }
    
    #pragma omp parallel private(sum) shared(total_sum)
    {

        if (omp_get_thread_num() == 0) {
            // printf("Node %d has available threads: %d\n", this_rank, omp_get_num_threads());
        }
        sum = 0;

        #pragma for
        for (int i = start; i < end; i++) {
            sum += my_func(i);
        }

        #pragma critical
        {
            total_sum += sum;
        }
    // for (int i = 0; i < N)

        if (omp_get_thread_num() == 0) 
            printf("Node %d computed sum: %lf with %d threads\n", this_rank, total_sum, omp_get_num_threads());
    }


    // Collect the results from each process

    MPI_Gather(&total_sum, 1, MPI_DOUBLE, received_sums->data, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (this_rank == 0) {
        printf("Rank 0 gathered matrix:\n");
        Matrix_print_d(received_sums);

        printf("sum: %lf\n", sum_d(received_sums));
    }

    MPI_Finalize();

    return 0;
}