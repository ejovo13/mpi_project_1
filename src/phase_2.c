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
int ntheta = -1;
char * binary_out = NULL;
char * precomputed_binary = NULL;

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

int main() {





    return 0;
}