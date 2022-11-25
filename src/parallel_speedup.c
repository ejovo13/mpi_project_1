/**========================================================================
 * ?                          parallel_speedup.c
 * @brief   : Time how effective the implementation of OpenMP parallelization
 * speeds up program execution
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-11-19
 *========================================================================**/

// We're gonna go ahead and collect speed up data for different model sizes, randing from 0 to 100 (to keep things fast)

// maybe 200 as well just to collect more data

#include "geodesy.h"

Matrix_d *time_model(int lmax);

#define MAX_THREADS 12

int main() {

    int threads[MAX_THREADS] = {0};
    for (int i = 0; i < MAX_THREADS; i++) { 
        threads[i] = i + 1;
    }

    Matrix_d *t_50 = time_model(50);
    Matrix_d *t_25 = time_model(25);
    Matrix_d *t_100 = time_model(100);
    Matrix_d *t_150 = time_model(150);

    // Let's write this data to a csv file
    FILE *csv_out = fopen("parallel_timing.csv", "w");

    fprintf(csv_out, "n_threads,t25,t50,t100,t150\n");
    // Let's iterate through the number of threads
    for (int i = 0; i < MAX_THREADS; i++) {
        fprintf(csv_out, "%d,%lf,%lf,%lf,%lf\n", threads[i], vecat_d(t_25, i), vecat_d(t_50, i), vecat_d(t_100, i), vecat_d(t_150, i));
    }

    // fprintf(csv_out, "n_threads,t25,t50\n");
    // // Let's iterate through the number of threads
    // for (int i = 0; i < MAX_THREADS; i++) {
    //     fprintf(csv_out, "%d,%lf,%lf\n", threads[i], vecat_d(t_25, i), vecat_d(t_50, i));
    // }


    return 0;
}

// Time how long it takes to compute a model for a given l, starting with 1 .. 24 threads
// We'll be using the small model as our data set for timing purposes.
// use the precomp for SUP_LMAX = 200
Matrix_d *time_model(int lmax) {

    const int SUP_LMAX = 200;

    int threads[MAX_THREADS] = {0};
    for (int i = 0; i < MAX_THREADS; i++) { 
        threads[i] = i + 1;
    }

    char plm_bin[100] = {0};
    sprintf(plm_bin, "ETOPO1_small_P%d.bin", SUP_LMAX);

    data_iso *data = get_data_small(true);
    // Precomputation data has already been computed for us
    Precomp *precomp = newPrecomp(0, lmax, SUP_LMAX, data, plm_bin);
    SphericalModel *model = newSphericalModel(lmax);


    // Time average execution for different threads
    int n_trials = 5;
    Matrix_d *avg_time = Matrix_new_d(1, MAX_THREADS);

    for (int i = 0; i < MAX_THREADS; i++) {

        double total_time = 0;

        for (int n = 0; n < n_trials; n++) {
            total_time += modelComputeCSlmPrecompOMP2Threads(model, data, precomp, threads[i]);
        }

        vecset_d(avg_time, i, total_time / n_trials);
    }

    // Return the matrix
    return avg_time;
}

