#include "mpi_util.h"

int compute_workload(int total_work, int world_size, int this_rank) {

    int workload = total_work / world_size;
    if (this_rank < (total_work % world_size)) workload ++;

    return workload;
}

Matrix_i *compute_workload_array(int total_work, const int world_size) {

    Matrix_i *workload_array = Vector_new_i(world_size);

    for (size_t i = 0; i < world_size; i++) {
        Vector_set_i(workload_array, i, compute_workload(total_work, world_size, i));
    }

    return workload_array;
}

Matrix_i *compute_displacements(int total_work, const int world_size) {

    // Create a new matrix length world_size
    Matrix_i *workload_array = compute_workload_array(total_work, world_size);
    Matrix_i *displacements = Vector_new_i(world_size);

    Vector_set_i(displacements, 0, 0);

    for (int i = 1; i < world_size; i++) {
        displacements->data[i] = displacements->data[i - 1] + (workload_array->data[i - 1]);
    }

    Matrix_free_i(workload_array);

    return displacements;

}

Matrix_i *compute_startend_array(int total_work, int world_size) {

    // Create a new matrix length world_size
    Matrix_i *workload_array = compute_workload_array(total_work, world_size);
    Matrix_i *startend = Vector_new_i(world_size + 1);

    // set the first element to -1
    Vector_set_i(startend, 0, -1);

    for (int i = 0; i < world_size; i++) {
        startend->data[i + 1] = startend->data[i] + workload_array->data[i];
    }

    Matrix_free_i(workload_array);

    return startend;
}