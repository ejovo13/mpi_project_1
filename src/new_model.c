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

int main(int argc, char **argv) {

    args_t *args = process_command_line_options(argc, argv, false);

    print_args(args);

    char plm_bin[100] = {0};
    char coeff_file_bin[100] = {0};

    sprintf(coeff_file_bin, "sph_%s_%d.bin", args->size_dataset, args->lmodel);
    sprintf(plm_bin, "ETOPO1_%s_P%d.bin", args->size_dataset, args->lbin);

    if (args->from) {
        assert(args->a != 0);
        printf("Computing model with degree %d from %d\n", args->lmodel, args->a);
    } else {
        printf("==============================\n");
        printf("Computing model of degree %d\n", args->lmodel);
        printf("==============================\n");
    }


    Precomp *precomp = newPrecomp(0, args->lmodel, args->lbin, args->data, plm_bin);
    SphericalModel *model = buildSphericalModel(args->data, args->lmodel, coeff_file_bin, args->recompute, args->from, args->a);

    printf("[main] Model built!\n");

    // output the model to a text file
    if (args->txt) {
        SphericalModelToTXT(model, args->size_dataset);
    }

    // compute the predicted values and store them in a binary file
    // fhat_<size_dataset>_<lmodel>.bin
    if (args->predict) {

        char binary_prediction_out[100] = {0};
        sprintf(binary_prediction_out, "fhat_%s_%d.bin", args->size_dataset, args->lmodel);
        Matrix_f *f_hat = NULL;
        Clock *clock = Clock_new();

        FILE *test_open = fopen(binary_prediction_out, "rb");
        if (test_open == NULL) {

            /**========================================================================
             *!                           Compute Prediction
            *========================================================================**/
            printf("[main] Computing altitude predictions...\n");
            Clock_tic(clock);
            f_hat = compute_prediction(model, precomp, args->data);
            Clock_toc(clock);
            printf("[main] Time to compute prediction: %lfs\n", elapsed_time(clock));
            save_prediction(f_hat, binary_prediction_out);
            Vector_print_head_f(f_hat, 10);

        } else {

            printf("[main] Predictions already computed\n");
            fclose(test_open);
        }

        head_data(args->data);

        // We can't reach this point before having the predictions already computed
        if (args->diff) {

            char diff_out[100] = {0};
            sprintf(diff_out, "diff_%s_%d.csv", args->size_dataset, args->lmodel);

            // test_open = fopen(diff_out, "r");
            // if (test_open == NULL) {

                Clock_tic(clock);
                f_hat = load_prediction(binary_prediction_out, args->data->N);
                Clock_toc(clock);
                printf("[main] Loaded predictions in %lfs\n", elapsed_time(clock));
                printf("[main] Computing differences\n");
                FILE *diff_csv = fopen(diff_out, "w");
                fprintf(diff_csv, "diff\n");
                Vector_print_head_f(f_hat, 10);

                for (int i = 0; i < args->data->N; i++) {
                    double res = f_hat->data[i] - (float) args->data->r[i];
                    printf("f_hat->data[i]: %f, data->r[i]: %u, residual: %f\n", f_hat->data[i], args->data->r[i], res);
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