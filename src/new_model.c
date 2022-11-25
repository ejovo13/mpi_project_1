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

    system("cat world_ascii.txt");
    printf("\n");

    args_t *args = process_command_line_options(argc, argv, false);

    if (args->print_args)
        print_args(args);

    char plm_bin[100] = {0};
    char coeff_file_bin[100] = {0};

    sprintf(coeff_file_bin, "sph_%s_%d.bin", args->size_dataset, args->lmodel);
    sprintf(plm_bin, "ETOPO1_%s_P%d.bin", args->size_dataset, args->lbin);

    if (args->from) {
        assert(args->a != 0);
        printf("Computing model with degree %d starting from %d\n", args->lmodel, args->a);
    } else {
        printf("==============================\n");
        printf("Computing model of degree %d\n", args->lmodel);
        printf("==============================\n");
    }

    char command[100] = {0};

    Precomp *precomp = newPrecomp(0, args->lmodel, args->lbin, args->data, plm_bin);
    SphericalModel *model = buildSphericalModel(
        args->data, 
        args->lmodel, 
        coeff_file_bin,
        args->recompute,
        args->from,
        args->a
    );

    printf("\n[main] Model built!\n");

    // output the model to a text file
    if (args->txt) {
        SphericalModelToTXT(model, args->size_dataset);
    }

    if (args->predict)
        predict_stuff(args, model, precomp);

}