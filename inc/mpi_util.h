/**========================================================================
 * ?                          mpi_util.h
 * @brief   : MPI utility functions to statically divide the workload
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-11-16
 *========================================================================**/
#include "ejovo.h"

int compute_workload(int total_work, int world_size, int this_rank);

Matrix_i *compute_workload_array(int total_work, const int world_size);

Matrix_i *compute_displacements(int total_work, const int world_size);

Matrix_i *compute_startend_array(int total_work, int world_size);