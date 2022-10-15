/**========================================================================
 * ?                          reduce.c
 * @brief   : Reduce the size of a ETOPO1 data set.
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-14
 *========================================================================**/
#include <stdio.h>
#include <stdlib.h>

#include "harmonics.h"

int main() {

    // reduce t p 
    // Reduce a data set and save it in the iso order

    const char filename[] = "ETOPO1_small.csv";
    // const char filename_out_tiny[] = "ETOP01_tiny.csv";
    // const char filename_out_[] = "ETOP01_tiny.csv";

    // load_data_iso(filename, )
    int initial_size = 64800;
    int p_reduction = 10;
    int t_reduction = 10;

    FILE *open = fopen(filename, "r");
    FILE *out_tiny  = fopen("ETOPO1_tiny.csv", "w");
    FILE *out_unit  = fopen("ETOPO1_unit.csv", "w");

    size_t buffsize = 200;
    // size_t character
    char *buffer = (char *) malloc(sizeof(*buffer) * buffsize);

    for (int i = 0; i < initial_size; i++) {
        getline(&buffer, &buffsize, open);
        if (i % 10 == 0 && (i / 180) % 10 == 0) {
            fprintf(out_tiny, "%s", buffer);
        }
        if (i % 30 == 0 && (i / 180) % 30 == 0) {
            fprintf(out_unit, "%s", buffer);
        }
    }

    printf("Running reduce\n");


    return 0;
}
