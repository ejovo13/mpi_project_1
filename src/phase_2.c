#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>

#include "geodesy.h"

/**========================================================================
 * ?                          phase_2.c
 * @brief   : Using precomputed Associated Legendre function data, pull data
 *            from a binary data set and compute the coefficients Clm and Slm
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-19
 *========================================================================**/

// I need the size of the data set 
int lmax = -1;
int lmodel = -1;
int ntheta = -1;
char * binary_out = NULL;
char * data_bin = NULL;
char * precomputed_binary = NULL;

void usage(char ** argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--ntheta n                      number of theta in discretization\n");
        printf("--lmax l                        degree of the stored binary file\n");
        printf("--lmodel L                      degree of the model to compute\n");
        printf("--data FILE                     name of the binary data file\n");
        printf("--in FILE                       name of the precomputed PLM values file\n");
        printf("--out FILE                      name of the output binary file\n");
        printf("\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[] = {
                {"ntheta", required_argument, NULL, 'n'},
                {"lmax", required_argument, NULL, 'l'},
                {"lmodel", required_argument, NULL, 'L'},
                {"data", required_argument, NULL, 'd'},
                {"out", required_argument, NULL, 'o'},
                {"in", required_argument, NULL, 'i'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'n':
                        ntheta = atoi(optarg);
                        break;
                case 'o':
                        binary_out = optarg;
                        break;
                case 'i':
                        precomputed_binary = optarg;
                        break;
                case 'd':
                        data_bin = optarg;
                        break;
                case 'L':
                        lmodel = atoll(optarg);
                        break;
                case 'l':
                        lmax = atoll(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        /* missing required args? */
        if (binary_out == NULL || precomputed_binary == NULL || data_bin == NULL || ntheta < 0 || lmax < 0 || lmodel < 0)
                usage(argv);
}

int main(int argc, char **argv) {

    process_command_line_options(argc, argv);

    data_iso *data = load_data_binary(data_bin, ntheta, ntheta * 2, true);

    head_data(data);

    Precomp *precomp = newPrecomp(0, lmodel, lmax, data, precomputed_binary);

    // Matrix_print_d(precomp->Plm_th);

    SphericalModel *model = newSphericalModel(lmodel);
    double time = modelComputeCSlmPrecomp(model, data, precomp);

    SphericalModelToTXT(model, "med");

    return 0;
}