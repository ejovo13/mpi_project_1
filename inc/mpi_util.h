/**========================================================================
 * ?                          mpi_util.h
 * @brief   : MPI utility functions to statically divide the workload
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-11-16
 *========================================================================**/
#include "ejovo.h"

// Run a block of code one time on the root process (rank 0). This macro relies on the 
// definition of <code>this_rank<\code>
#define MPI_ONCE(block) if (this_rank == 0) {\
    block \
    }

// Run a block of code in order for i = [0, 1, ..., world_size - 1]
#define MPI_ORDERED(block) for (int i = 0; i < world_size; i++) {\
    if (this_rank == i) { \
        block \
    } \
    MPI_Barrier(MPI_COMM_WORLD); \
}

int compute_workload(int total_work, int world_size, int this_rank);

Matrix_i *compute_workload_array(int total_work, const int world_size);

Matrix_i *compute_displacements(int total_work, const int world_size);

Matrix_i *compute_startend_array(int total_work, int world_size);