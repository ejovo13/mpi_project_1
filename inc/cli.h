#ifndef GEODESY_CLI_H
#define GEODESY_CLI_H

/**========================================================================
 * ?                          cli.h
 * @brief   : Manage CLI parsing for different models
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-11-25
 *========================================================================**/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>

#include "data.h"
#include "phase_1.h"
#include "phase_2.h"


/**========================================================================
 *!                           Sequential
 *========================================================================**/
typedef struct {

    int argc;
    char ** argv;
    int lmodel;
    int lbin;
    char *size_dataset;
    data_iso *data;

    bool from;
    int a;

    bool txt;
    bool predict;
    bool diff;
    bool recompute;
    bool print_args;

    const char *plm_bin;
    const char *coeff_file_bin;

} args_t;

// bool to string
static inline const char *btos(bool b) {
    if (b) return "true";
    else return "false";
}

void print_args(const args_t *args);

args_t *new_args(int argc, char **argv);

void usage(char ** argv);

args_t *process_command_line_options(int argc, char ** argv, bool __log);

#endif // GEODESY_CLI_