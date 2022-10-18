#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>

#include "geodesy.h"

int lmax = -1;
int npoint = -1;
char * data_filename;
char * model_filename;

void usage(char ** argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--data FILENAME              input file containing experimental data points\n");
        printf("--model FILENAME             output file containing the model\n");
        printf("--npoint N                   number of points to read\n");
        printf("--lmax N                     order of the model\n");
        printf("\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[5] = {
                {"data", required_argument, NULL, 'd'},
                {"npoint", required_argument, NULL, 'n'},
                {"lmax", required_argument, NULL, 'l'},
                {"model", required_argument, NULL, 'm'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'd':
                        data_filename = optarg;
                        break;
                case 'm':
                        model_filename = optarg;
                        break;
                case 'n':
                        npoint = atoi(optarg);
                        break;
                case 'l':
                        lmax = atoll(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        /* missing required args? */
        if (data_filename == NULL || model_filename == NULL || lmax < 0 || npoint <= 0)
                usage(argv);
}

int main(int argc, char ** argv)
{
	process_command_line_options(argc, argv);
    old_model(lmax, npoint, data_filename, model_filename);
}