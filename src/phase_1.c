#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>

/**========================================================================
 * ?                          phase_1.c
 * @brief   : Precompute Associated Legendre values Plm(cos\theta) for all
 *            values of theta for a dataset given ntheta
 * 
 * @example: 
 *            ./phase1 --ntheta 180 --lmax 500 --out ETOPO1_small_P500.bin         
 *  
 * 
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-19
 *========================================================================**/

#include "geodesy.h"

int lmax = -1;
int ntheta = -1;
char * binary_out = NULL;

void usage(char ** argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--ntheta N                      number of theta in discretization\n");
        printf("--lmax L                        degree of the stored binary file\n");
        printf("--out FILE                      name of the output binary file\n");
        printf("\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[] = {
                {"ntheta", required_argument, NULL, 'n'},
                {"lmax", required_argument, NULL, 'l'},
                {"out", required_argument, NULL, 'o'},
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
                case 'l':
                        lmax = atoll(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        /* missing required args? */
        if (binary_out == NULL || ntheta < 0 || lmax < 0)
                usage(argv);
}


int main(int argc, char **argv) {

    process_command_line_options(argc, argv);
    // Get the number of theta, the name of the model, and the max l value to compute a binary file
    // storing the contents of P_lm_th

    const double d_theta = PI / ntheta;

    // Let's create a theta matrix 
    Matrix_d *theta = Matrix_new_d(1, ntheta);

    for (int i = 0; i < ntheta; i++) {
        theta->data[i] = d_theta * i;
    }


    write_binary_plm(lmax, theta, binary_out); // Runs in about 30s for the Ultra data set and L = 788


    return 0;
}