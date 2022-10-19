#include "geodesy.h"

// In phase_1 and phase_2 we created a binary file that has precomputations in it, and we also computed a new 
// spherical model.

// In this phase we want to load in the model and compare the values with a particulary data set.

int main() {

    const char *binary_plm  = "ETOPO1_small_P1000.bin";
    const char *binary_cslm = "sph_100_small.bin"; 
    const char *binary_prediction_out = "fhat_100_small.bin";

    data_iso *data = get_data_small();

    const int L0 = 0;
    const int LF = 100;
    const int LMAX_BIN = 1000;

    Clock *clock = Clock_new();
    Precomp *precomp = newPrecomp(L0, LF, LMAX_BIN, data, binary_plm); 
    Matrix_f *f_hat = NULL;

    /**========================================================================
     *!                           Compute model
     *========================================================================**/
    // SphericalModel *model = newSphericalModel(LF, data, precomp);
    // double time = modelComputeCSlmPrecomp(model, data, precomp);
    // printf("Time to compute new model at LF %d: %lfs\n", LF, time);

    /**========================================================================
     *!                        or  Load Model
     *========================================================================**/
    Clock_tic(clock);
    SphericalModel *model = loadSphericalModel(binary_cslm, LF);
    Clock_toc(clock);
    printf("Loaded model in %lfs\n", elapsed_time(clock));

    /**========================================================================
     *!                           Compute Prediction
     *========================================================================**/
    // Clock_tic(clock);
    // f_hat = compute_prediction(model, precomp, data);
    // Clock_toc(clock);
    // printf("Time to compute prediction: %lfs\n", elapsed_time(clock));
    // save_prediction(f_hat, binary_prediction_out);
    // Matrix_free_f(f_hat);


    /**========================================================================
     *!                        or Load Prediction
     *========================================================================**/
    Clock_tic(clock);
    f_hat = load_prediction(binary_prediction_out, data->N);
    Clock_toc(clock);
    printf("Loaded prediction in %lfs\n", elapsed_time(clock));

    Vector_print_head_f(f_hat, 10);

    for (int i = 0; i < 20; i++) {
        printf("%f ", f_hat->data[i]);
    }
    printf("\n");

    // Now go through and compute the mse.
    double abs_error = 0;
    double sq_error  = 0;
    double err       = 0;
    // for (int i = 0; i < data->N; i++) {
    for (int i = 0; i < data->N; i++) {
        // printf("i: %d (%d - %f)\n", i, data->r[i], f_hat->data[i]);
        err = ((float) data->r[i]) - f_hat->data[i];
        abs_error += fabs(err);
        sq_error  += err * err;
    }

    printf("abs_error:    %lf\n", abs_error);
    printf("sq_error:     %lf\n", sq_error);
    printf("MSE:          %lf\n", sq_error / data->N);
    printf("Average err:  %lf\n", abs_error / data->N);

    return 0;
}