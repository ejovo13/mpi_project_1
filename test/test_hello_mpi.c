/**========================================================================
 * ?                          test_hellp_mpi.c
 * @brief   : Test mpi with address sanitizing to see if there are errors or not
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-11-26
 *========================================================================**/
#include <stdio.h>
#include <mpi.h>

int main(int argc, char ** argv) {

    int world_size, this_rank; 


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

    printf("Sup yo");


    // __attribute__((no_sanitize_address)) MPI_Finalize();
    MPI_Finalize();
    return 0;
}