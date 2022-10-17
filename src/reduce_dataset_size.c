/**========================================================================
 * ?                          reduce_dataset_size.c
 * @brief   : Reduce the size of a data set by sotring the z values
 *            as int16_t in a binary file
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-17
 *========================================================================**/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <stdint.h>

#include <geodesy.h>

int npoint = -1;
char * out = NULL;
char * in = NULL;

void usage(char ** argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--in FILENAME               input file containing experimental data points\n");
        printf("--out FILENAME              output file containing binary z values\n");
        printf("--npoint N                      number of theta in discretization\n");
        printf("\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[] = {
                {"npoint", required_argument, NULL, 'n'},
                {"out", required_argument, NULL, 'o'},
                {"in", required_argument, NULL, 'i'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'n':
                        npoint = atoi(optarg);
                        break;
                case 'i':
                        in = optarg;
                        break;
                case 'o':
                        out = optarg;
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        /* missing required args? */
        if (in == NULL || out == NULL || npoint < 0)
                usage(argv);
}

int main(int argc, char ** argv) {

    process_command_line_options(argc, argv);

    printf("Processing %d points of data\n", npoint);

    reduce_csv_to_binary(npoint, in, out);

    return 0;
}

