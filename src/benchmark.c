/**========================================================================
 * ?                          benchmark.c
 * @brief   : Gather timing data on the different models
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-18
 *========================================================================**/

#include "geodesy.h"

int main() {

    // It's ugly, but some of these file paths are going to be temporarily hardcoded. I will work on making this
    // more portable at a later date.

    // Run timing results on the small data set to start

    // We are assuming that this file is being executed in the project_home/build/src directory.
    const char *data_file = "../../csv/ETOPO1_small.csv";
    const char *model_file = "tmp";
    const int npoint = 64800;

    // Start with timing results for old model

    // Let's start of by only running each test 10 times, from l = 0 to 10
    const int L0 = 0;
    // const int LMAX = 5;  // like 5 seconds with nruns = 5
    // const int LMAX = 10; // 49s with nruns = 5
    const int LMAX = 20;  // 830 seconds with nruns = 5. Damn
    const int nruns = 5;

    // Matrix_d *tmp_row = Matrix_new_d(1, nruns);
    Matrix_d *old_benchmark = Matrix_new_d(LMAX + 1, 1);
    Matrix_d *new_benchmark = Matrix_new_d(LMAX + 1, 1);

    // Time old model
    for (int i = L0; i <= LMAX; i++) {
        double sum = 0;
        for (int n = 0; n < nruns; n++) {
            sum += old_model(i, npoint, data_file, model_file);            
        }
        old_benchmark->data[i] = sum / nruns;
    }

    // Time new model
    for (int i = L0; i <= LMAX; i++) {
        double sum = 0;
        for (int n = 0; n < nruns; n++) {
            sum += time_new_model(i, npoint, data_file);            
            // sum += time_new_model(i, npoint, data_file, model_file);            
        }
        new_benchmark->data[i] = sum / nruns;
    }

    // now I need to actually output this data to a csv.

    FILE *timing_csv = fopen("timing.csv", "w");

    // Now loop through the rows
    fprintf(timing_csv, "l,old,new\n");

    for (int i = L0; i <= LMAX; i++) {
        fprintf(timing_csv, "%d,%lf,%lf\n", i, vecat_d(old_benchmark, i), vecat_d(new_benchmark, i));
    }

    fclose(timing_csv);

    return 0;
}