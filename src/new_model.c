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
SphericalModel *build_model(const args_t *args);

int main(int argc, char **argv) {

    args_t *args = process_command_line_options(argc, argv, false, 0);

    if (args->ascii) {
        system("cat world_ascii.txt");
        printf("\n");
    }

    if (args->help) {
        usage(argv);
        exit(0);
    }

    if (args->print_args)
        print_args(args);

    if (args->from) {
        assert(args->a != 0);
        printf("Computing model with degree %d starting from %d\n", args->lmodel, args->a);
    } else {
        printf("==============================\n");
        printf("Computing model of degree %d\n", args->lmodel);
        printf("==============================\n");
    }

    SphericalModel *model = build_model(args); 

    printf("[main] Model succesfully built!\n");

    // output the model to a text file
    if (args->txt) {
        SphericalModelToTXT(model, args->size_dataset);
    }

    if (args->predict) {
        Precomp *precomp = newPrecomp(0, args->lmodel, args->lbin, args->data, args->plm_bin);
        predict_stuff(args, model, precomp);
        freePrecomp(precomp);
    }

    freeSphericalModel(model);
    freeArgs(args);

}

SphericalModel *build_model(const args_t *args) {
    return buildSphericalModel(
        args->data, 
        args->lmodel, 
        args->coeff_file_bin,
        args->recompute,
        args->from,
        args->a
    );
}