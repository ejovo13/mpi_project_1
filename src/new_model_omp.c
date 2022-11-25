#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>

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
int lmax = -1;
int lmodel = -1;
char * size_dataset = NULL;
data_iso *data = NULL;

bool txt = false;
bool predict = false;
bool diff = false;

void usage(char ** argv)
{
    printf("%s [OPTIONS]\n\n", argv[0]);
    printf("Options:\n");
    printf("--[small | med | hi | ultra]    dataset to model\n");
    printf("--lmodel l                      degree of the model to compute\n");
    printf("--lmax L                        degree of the stored binary file\n");
    printf("[--txt]                         output the model as a text file\n");
    printf("[--predict]                     predict the altitude values (f_hat)\n");
    printf("[--diff]                        predict the altitude values (f_hat) and compute\n");
    printf("                                the difference between the model and predicted value\n");
    printf("\n");
    exit(0);
}

void process_command_line_options(int argc, char ** argv)
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
        {NULL, 0, NULL, 0}
    };
    char ch;
    while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch (ch) {

        case 'L':
            lmodel = atoll(optarg);
            break;
        case 'l':
            lmax = atoll(optarg);
            break;
        case 's':
            size_dataset = "small";
            data = get_data_small(true);
            break;
        case 'm':
            size_dataset = "med";
            data = get_data_med(true);
            break;
        case 'h':
            size_dataset = "hi";
            data = get_data_hi(true);
            break;
        case 'u':
            size_dataset = "ultra";
            data = get_data_ultra(true);
            break;
        case 't':
            txt = true;
            break;
        case 'p':
            predict = true;
            break;
        case 'd':
            diff = true;
            predict = true;
            break;
        default:
            errx(1, "Unknown option\n");
        }
    }
    /* missing required args? */
    if (size_dataset == NULL || data == NULL || lmax < 0 || lmodel < 0)
        usage(argv);
}

int main(int argc, char **argv) {

    process_command_line_options(argc, argv);

    char plm_bin[100] = {0};
    char coeff_file_bin[100] = {0};
    sprintf(coeff_file_bin, "sph_%s_%d.bin", size_dataset, lmodel);
    sprintf(plm_bin, "ETOPO1_%s_P%d.bin", size_dataset, lmax);

    Precomp *precomp = newPrecomp(0, lmodel, lmax, data, plm_bin);
    SphericalModel *model = NULL;

    FILE *test_open = fopen(coeff_file_bin, "rb");
    if (test_open == NULL) {

        printf("[main] %s not found, computing Clm and Slm coefficients\n", coeff_file_bin);
        model = newSphericalModel(lmodel);   
        double time = modelComputeCSlmPrecompOMP(model, data, precomp);
        printf("[main] Computed coefficients for L = %d in %lfs\n", lmodel, time);
        SphericalModelToBIN(model, size_dataset);

    } else {

        // Model file coeff_file_bin already exists, so load from it
        fclose(test_open);
        model = loadSphericalModel(coeff_file_bin, lmodel);

    }

    // output the model to a text file
    if (txt) {
        SphericalModelToTXT(model, size_dataset);
    }

    // compute the predicted values and store them in a binary file
    // fhat_<size_dataset>_<lmodel>.bin
    if (predict) {

        char binary_prediction_out[100] = {0};
        sprintf(binary_prediction_out, "fhat_%s_%d.bin", size_dataset, lmodel);
        Matrix_f *f_hat = NULL;
        Clock *clock = Clock_new();

        FILE *test_open = fopen(binary_prediction_out, "rb");
        if (test_open == NULL) {

            /**========================================================================
             *!                           Compute Prediction
                *========================================================================**/
            printf("[main] Computing altitude predictions...\n");
            Clock_tic(clock);
            f_hat = compute_prediction(model, precomp, data);
            Clock_toc(clock);
            printf("[main] Time to compute prediction: %lfs\n", elapsed_time(clock));
            save_prediction(f_hat, binary_prediction_out);
            Vector_print_head_f(f_hat, 10);

        } else {

            printf("[main] Predictions already computed\n");
            fclose(test_open);
        }

        head_data(data);

        // We can't reach this point before having the predictions already computed
        if (diff) {

            char diff_out[100] = {0};
            sprintf(diff_out, "diff_%s_%d.csv", size_dataset, lmodel);

            // test_open = fopen(diff_out, "r");
            // if (test_open == NULL) {

                Clock_tic(clock);
                f_hat = load_prediction(binary_prediction_out, data->N);
                Clock_toc(clock);
                printf("[main] Loaded predictions in %lfs\n", elapsed_time(clock));
                printf("[main] Computing differences\n");
                FILE *diff_csv = fopen(diff_out, "w");
                fprintf(diff_csv, "diff\n");
                Vector_print_head_f(f_hat, 10);

                for (int i = 0; i < data->N; i++) {
                    double res = f_hat->data[i] - (float) data->r[i];
                    // printf("f_hat->data[i]: %f, data->r[i]: %u, residual: %d\n", f_hat->data[i], data->r[i], res);
                    fprintf(diff_csv, "%f\n", res);
                }

                fclose(diff_csv);

            // } else {

                // printf("[main] Differences already computed, stored in %s\n", diff_out);
                // fclose(test_open);
            // } // if (test_open == NULL)
        } // if (diff)
    } // if (predict)
}