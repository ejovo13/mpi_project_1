#include "geodesy.h"

/**========================================================================
 * ?                          new_model.c
 * @brief   : Compute a Spherical Model of the earth when provided an L-value
 *            and the size of the data set
 * 
 * @usage   : ./new_model <SIZE_DATASET> --lmodel <L_MODEL> [--lmax <MAX_L>]  
 * @example : ./new_model --small --lmodel 20 --lmax 1000
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-28
 *========================================================================**/
// Arguments

SphericalModel *build_model_mpi(const args_t *args, int world_size, int this_rank);

int main(int argc, char **argv) {

    int this_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

    MPI_Barrier(MPI_COMM_WORLD);

    bool __log = false;


    if (this_rank == 0) {
        printf("Running mpi version with %d processes\n", world_size);
        __log = true;
    }

    args_t *args = process_command_line_options(argc, argv, __log, this_rank);
    args->rank = this_rank;
    args->world_size = world_size;

    if (args->ascii) {
        if (this_rank == 0) {
            system("cat world_ascii.txt");
            printf("\n");
        }
    }

    if (args->help) {
        if (this_rank == 0) 
            usage(argv);
        MPI_Finalize();
        exit(0);
    }

    if (this_rank == 0 && args->print_args) {
        print_args(args);
    }


    MPI_ONCE(
        if (args->from) {
            assert(args->a != 0);
            printf("Computing model with degree %d starting from %d\n", args->lmodel, args->a);
        } else {
            printf("==============================\n");
            printf("Computing model of degree %d\n", args->lmodel);
            printf("==============================\n");
        }
    )

    char plm_bin[100] = {0};
    char coeff_file_bin[100] = {0};
    
    sprintf(coeff_file_bin, "sph_%s_%d.bin", args->size_dataset, args->lmodel);
    sprintf(plm_bin, "ETOPO1_%s_P%d.bin", args->size_dataset, args->lbin);

    SphericalModel *model = build_model_mpi(args, world_size, this_rank);

    assert(model);
    assert(args->data);

    if (this_rank == 0) printf("[main] MPI model built\n");

    // output the model to a text file
    if (args->txt) {
        if (this_rank == 0)
            SphericalModelToTXT(model, args->size_dataset);
    }

    if (args->predict) {

        Precomp *precomp = newPrecomp(0, args->lmodel, args->lbin, args->data, args->plm_bin);
        predict_stuff_mpi(args, model, precomp);
        freePrecomp(precomp);
    }

    if (this_rank == 0) {
        printf("Reached finalize\n");
    }

    freeSphericalModel(model);
    freeArgs(args);

    MPI_Finalize();
    return 0;
}

SphericalModel *build_model_mpi(const args_t *args, int world_size, int this_rank) {
    return buildSphericalModelMPI(args->data, args->lmodel, args->lbin, args->coeff_file_bin, args->recompute, args->from, args->a, world_size, this_rank);
}

