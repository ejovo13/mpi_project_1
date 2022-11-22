#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <stdio.h>
#include <mpi.h>

#include "geodesy.h"
#include "ejovo.h"
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
typedef struct mpi_args {

    int lmax;
    int lmodel;
    char * size_dataset;
    data_iso *data;

    bool txt;
    bool predict;
    bool diff;
    bool recompute;

} mpi_args;

void init_mpi_args(mpi_args *args) {
    args->lmax = -1;
    args->lmodel = -1;
    args->size_dataset = NULL;
    args->data = NULL;
    args->txt = false;
    args->predict = false;
    args->diff = false;
    args->recompute = false;
}

void print_mpi_args(mpi_args *args) {
    printf("args = mpi_args {.lmax = %d, .lmodel = %d, .size_dataset = %s, .txt = %d, .predict = %d, diff = %d, r = %d\n", 
        args->lmax,
        args->lmodel,
        args->size_dataset,
        args->txt,
        args->predict,
        args->diff,
        args->recompute);
}

void usage_mpi(char ** argv, int this_rank)
{

    if(this_rank == 0) {

        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--[small | med | hi | ultra]    dataset to model\n");
        printf("--lmodel l                      degree of the model to compute\n");
        printf("--lmax L                        degree of the stored binary file\n");
        printf("[--txt]                         output the model as a text file\n");
        printf("[--predict]                     predict the altitude values (f_hat)\n");
        printf("[--diff]                        predict the altitude values (f_hat) and compute\n");
        printf("                                the difference between the model and predicted value\n");
        printf("[--recompute]                   compute the model coefficients no matter what files\n");
        printf("                                are present\n");
        printf("\n");

    }

    MPI_Finalize();
    exit(0);
}

void process_command_line_options(int argc, char ** argv, int this_rank, mpi_args *args)
{
    struct option longopts[] = {
        {"lmax", required_argument, NULL, 'l'},
        {"lmodel", required_argument, NULL, 'L'},
        {"small", no_argument, NULL, 's'},
        {"med", no_argument, NULL, 'm'},
        {"hi", no_argument, NULL, 'h'},
        {"ultra", no_argument, NULL, 'u'},
        {"txt", no_argument, NULL, 't'},
        {"predict", no_argument, NULL, 'p'},
        {"diff", no_argument, NULL, 'd'},
        {"recompute", no_argument, NULL, 'r'},

        {NULL, 0, NULL, 0}
    };

    init_mpi_args(args);

    char ch;

    while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch (ch) {

        case 'L':
            args->lmodel = atoll(optarg);
            break;
        case 'l':
            args->lmax = atoll(optarg);
            break;
        case 's':
            args->size_dataset = "small";
            args->data = get_data_small();
            break;
        case 'm':
            args->size_dataset = "med";
            args->data = get_data_med();
            break;
        case 'h':
            args->size_dataset = "hi";
            args->data = get_data_hi();
            break;
        case 'u':
            args->size_dataset = "ultra";
            args->data = get_data_ultra();
            break;
        case 't':
            args->txt = true;
            break;
        case 'p':
            args->predict = true;
            // printf("[USAGE]: predict set to %d\n", args->predict);
            break;
        case 'd':
            args->diff = true;
            args->predict = true;
            break;
        case 'r':
            args->recompute = true;
            break;
        default:
            errx(1, "Unknown option\n");
        }
    }
    /* missing required args? */
    if (args->size_dataset == NULL || args->data == NULL || args->lmax < 0 || args->lmodel < 0) 
        usage_mpi(argv, this_rank);

    if (args->size_dataset == NULL) {
        printf("Missing dataset size [--small | --med | --hi | --ultra] \n");
        usage_mpi(argv, this_rank);
    }

        
}



int main(int argc, char **argv) {

    int this_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

    if (this_rank == 0) {
        printf("Running mpi version with %d processes\n", world_size);
    }

    mpi_args args; 
    process_command_line_options(argc, argv, this_rank, &args);

    if (this_rank == 0) {
        print_mpi_args(&args);
    }


    char plm_bin[100] = {0};
    char coeff_file_bin[100] = {0};
    sprintf(coeff_file_bin, "sph_%s_%d.bin", args.size_dataset, args.lmodel);
    sprintf(plm_bin, "ETOPO1_%s_P%d.bin", args.size_dataset, args.lmax);

    Precomp *precomp = NULL;
    SphericalModel *model = NULL;

    MPI_ORDERED (
        precomp = newPrecomp(0, args.lmodel, args.lmax, args.data, plm_bin);
    )
    
    FILE *test_open = fopen(coeff_file_bin, "rb");
    if (test_open == NULL || args.recompute) {

        MPI_ONCE (
            printf("[main] %s not found, computing Clm and Slm coefficients\n", coeff_file_bin);
        )
        model = newSphericalModel(args.lmodel, args.data, precomp);    

        MPI_Barrier(MPI_COMM_WORLD);
        
        double time = modelComputeCSlmPrecompMPI(model, args.data, precomp, world_size, this_rank);

        // only write the model to a binary file if we are the 0th rank
        MPI_ONCE (
            printf("[main] Computed coefficients for L = %d in %lfs\n", args.lmodel, time);
            SphericalModelToBIN(model, args.size_dataset);
        )

    } else {

        // Model file coeff_file_bin already exists, so load from it
        fclose(test_open);
        model = loadSphericalModel(coeff_file_bin, args.lmodel);

    }

    // output the model to a text file
    if (args.txt) {
        if (this_rank == 0)
            SphericalModelToTXT(model, args.size_dataset);
    }

    // compute the predicted values and store them in a binary file
    // fhat_<size_dataset>_<lmodel>.bin
    if (args.predict) {

        char binary_prediction_out[100] = {0};
        sprintf(binary_prediction_out, "fhat_%s_%d.bin", args.size_dataset, args.lmodel);
        Matrix_f *f_hat = NULL;
        Clock *clock = Clock_new();

        FILE *test_open = fopen(binary_prediction_out, "rb");
        if (test_open == NULL) {

            /**========================================================================
             *!                           Compute Prediction
                *========================================================================**/
            printf("[main] Computing altitude predictions...\n");
            Clock_tic(clock);
            f_hat = compute_prediction_omp(model, precomp, args.data);
            Clock_toc(clock);
            printf("[main] Time to compute prediction: %lfs\n", elapsed_time(clock));
            save_prediction(f_hat, binary_prediction_out);
            Vector_print_head_f(f_hat, 10);

        } else {

            printf("[main] Predictions already computed\n");
            fclose(test_open);
        }

        head_data(args.data);

        // We can't reach this point before having the predictions already computed
        if (args.diff) {

            char diff_out[100] = {0};
            sprintf(diff_out, "diff_%s_%d.csv", args.size_dataset, args.lmodel);

            // test_open = fopen(diff_out, "r");
            // if (test_open == NULL) {

                Clock_tic(clock);
                f_hat = load_prediction(binary_prediction_out, args.data->N);
                Clock_toc(clock);
                printf("[main] Loaded predictions in %lfs\n", elapsed_time(clock));
                printf("[main] Computing differences\n");
                FILE *diff_csv = fopen(diff_out, "w");
                fprintf(diff_csv, "diff\n");
                printf("Loaded as f_hat: \n");

                Matrix_d *true_value = Matrix_new_d(1, args.data->N);

                for (int i = 0; i < args.data->N; i++) {
                    vecset_d(true_value, i, args.data->r[i]);
                }

                Matrix_d *diff_vec = Matrix_new_d(1, args.data->N);

                for (int i = 0; i < args.data->N; i++) {
                    float delta = f_hat->data[i] - (float) args.data->r[i];
                    vecset_d(diff_vec, i, delta);
                    // printf("f_hat->data[%d]: %f, data->r[i]: %u, diff: %lf\n", i, f_hat->data[i], data->r[i], diff);
                    fprintf(diff_csv, "%lf\n", delta);
                }

                Vector_print_head_f(f_hat, 10);
                Vector_print_head_d(true_value, 10);
                Vector_print_head_d(diff_vec, 10);

                fclose(diff_csv);

            // } else {

                // printf("[main] Differences already computed, stored in %s\n", diff_out);
                // fclose(test_open);
            // } // if (test_open == NULL)
        } // if (diff)
    } // if (predict)

    MPI_Finalize();
    return 0;
}