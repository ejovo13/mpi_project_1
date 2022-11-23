#include "geodesy.h"

// In phase_1 and phase_2 we created a binary file that has precomputations in it, and we also computed a new 
// spherical model.


typedef struct {

    double avg_error;
    double max_error;
    double min_error;
    double mse;
    double total_error;

} ErrorReport;

// typedef struct {

//     int sample_size;
//     ErrorReport; 

// } StochasticErrorReport;

// assuming the file is opened, add a line
void add_line_error_report(int l, ErrorReport report, FILE *file) {
    fprintf(file, "%d,%lf,%lf,%lf,%lf,%lf\n", l, report.avg_error, report.max_error, report.min_error, report.mse, report.total_error);
}



// In this phase we want to load in the model and compare the values with a particulary data set.
// Return a Matrix_d in the form [avg_error, max_error, min_error, mse, total_error]
ErrorReport t_validate_predictions_small(int l_model);

// Simply load the f_hat of a model of degree `l_model`, compute the difference between
// the prediction and the true function value, then output the difference to a csv file
// for statistical analysis in R
void t_write_predictions_small(int l_model); 





int main() {

    const int indices_array[] = {5, 10, 20, 30, 50, 100, 150, 175, 190, 195, 196, 197, 198, 199, 200, 201, 250, 300, 400, 500};

    // // Let the compiler help you out with counting this one
    const int n_indices = 20;

    Matrix_i *ind = Matrix_from_i(indices_array, 1, n_indices);
    // ErrorReport err_report;

    // FILE *err_report_out = fopen("error_report.csv", "w");

    // fprintf(err_report_out, "l,avg_error,max_error,min_error,mse,total_error\n");
    // for (int i = 0; i < n_indices; i++) {
    //     err_report = t_validate_predictions_small(ind->data[i]);
    //     add_line_error_report(ind->data[i], err_report, err_report_out);
    // }

    // fclose(err_report_out);


    // Now go ahead and write prediction files in csv format

    for (int i = 0; i < n_indices; i++) {
        t_write_predictions_small(ind->data[i]);
    }

    return 0;
}

ErrorReport t_validate_predictions_small(int l_model) {

    // Always use the P1000 data set and small resolution
    const char *binary_plm  = "ETOPO1_small_P1000.bin";
    // const char *binary_cslm = "sph_100_small.bin"; 
    // const char *binary_prediction_out = "fhat_100_small.bin";

    char binary_cslm[100] = {0};
    char binary_prediction_out[100] = {0};

    sprintf(binary_cslm, "sph_%d_small.bin", l_model);
    sprintf(binary_prediction_out, "fhat_%d_small.bin", l_model);

    data_iso *data = get_data_small(true);

    const int L0 = 0;
    const int LF = l_model;
    const int LMAX_BIN = 1000;

    Clock *clock = Clock_new();
    Precomp *precomp = newPrecomp(L0, LF, LMAX_BIN, data, binary_plm); 
    Matrix_f *f_hat = NULL;

    // The model is supposedly being computed in phase_1. This code is redundant
    // TODO remove compute model section
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


    FILE *test_open = fopen(binary_prediction_out, "rb");
    if (test_open == NULL) {

        /**========================================================================
         *!                           Compute Prediction
            *========================================================================**/
        Clock_tic(clock);
        f_hat = compute_prediction(model, precomp, data);
        Clock_toc(clock);
        printf("Time to compute prediction: %lfs\n", elapsed_time(clock));
        save_prediction(f_hat, binary_prediction_out);

    } else {

        fclose(test_open);
        /**========================================================================
         *!                        or Load Prediction
        *========================================================================**/
        Clock_tic(clock);
        f_hat = load_prediction(binary_prediction_out, data->N);
        Clock_toc(clock);
        printf("Loaded prediction in %lfs\n", elapsed_time(clock));

        Vector_print_head_f(f_hat, 5);

    }

    // Now go through and compute the mse.
    double abs_error = 0;
    double sq_error  = 0;
    double err       = 0;
    double min_error = fabs((float) data->r[0] - f_hat->data[0]);
    double max_error = fabs((float) data->r[0] - f_hat->data[0]);
    // for (int i = 0; i < data->N; i++) {
    for (int i = 0; i < data->N; i++) {
        // printf("i: %d (%d - %f)\n", i, data->r[i], f_hat->data[i]);
        err = ((float) data->r[i]) - f_hat->data[i];
        abs_error += fabs(err);
        sq_error  += err * err;
        if (fabs(err) < min_error) min_error = fabs(err);
        if (fabs(err) > max_error) max_error = fabs(err);
    }

    printf("abs_error:    %lf\n", abs_error);
    printf("sq_error:     %lf\n", sq_error);
    printf("MSE:          %lf\n", sq_error / data->N);
    printf("Average err:  %lf\n", abs_error / data->N);

    ErrorReport report;

    report.avg_error = abs_error / data->N;
    report.max_error = max_error;
    report.min_error = min_error;
    report.mse = sq_error / data->N;
    report.total_error = abs_error;

    freePrecomp(precomp);
    freeSphericalModel(model);
    Matrix_free_f(f_hat);

    return report;
}

void t_write_predictions_small(int l_model) {

    // Load data
    const data_iso *data = get_data_small(true);

    // Load f_hat
    char binary_prediction_in[100] = {0};
    char diff_out[100] = {0};

    sprintf(binary_prediction_in, "fhat_%d_small.bin", l_model);
    sprintf(diff_out, "diff_%d_small.csv", l_model);
    Matrix_f *f_hat = load_prediction(binary_prediction_in, data->N);

    // Compute the difference and store it in a matrix
    float diff = 0;

    FILE *diff_csv = fopen(diff_out, "w");
    fprintf(diff_csv, "diff\n");

    for (int i = 0; i < data->N; i++) {
        diff = f_hat->data[i] - (float) data->r[i];
        fprintf(diff_csv, "%f\n", diff);
    }

    fclose(diff_csv);

}


// Let's leave this alone for now. We are going to perform statistic analysis in R.
// StochasticErrorReport t_stochastic_predictions_small(int l_model, int n) {

//     // Always use the P1000 data set and small resolution
//     const char *binary_plm  = "ETOPO1_small_P1000.bin";
//     // const char *binary_cslm = "sph_100_small.bin"; 
//     // const char *binary_prediction_out = "fhat_100_small.bin";

//     char binary_cslm[100] = {0};
//     char binary_prediction_out[100] = {0};

//     sprintf(binary_cslm, "sph_%d_small.bin", l_model);
//     sprintf(binary_prediction_out, "fhat_%d_small.bin", l_model);

//     data_iso *data = get_data_small();

//     const int L0 = 0;
//     const int LF = l_model;
//     const int LMAX_BIN = 1000;

//     Clock *clock = Clock_new();
//     Precomp *precomp = newPrecomp(L0, LF, LMAX_BIN, data, binary_plm); 
//     Matrix_f *f_hat = NULL;

//     /**========================================================================
//      *!                   Load computed model, produce estimations
//      *========================================================================**/
//     Clock_tic(clock);
//     SphericalModel *model = loadSphericalModel(binary_cslm, LF);
//     Clock_toc(clock);
//     printf("Loaded model in %lfs\n", elapsed_time(clock));

//     //! Ok, so our predictions have actually already been computed, and this
//     //! test file is just for seeing how the statistics behave
//     Clock_tic(clock);
//     f_hat = compute_prediction(model, precomp, data);
//     Clock_toc(clock);
//     printf("Time to compute prediction: %lfs\n", elapsed_time(clock));
//     save_prediction(f_hat, binary_prediction_out);




//     // Now go through and compute the mse.
//     double abs_error = 0;
//     double sq_error  = 0;
//     double err       = 0;
//     double min_error = fabs((float) data->r[0] - f_hat->data[0]);
//     double max_error = fabs((float) data->r[0] - f_hat->data[0]);
//     // for (int i = 0; i < data->N; i++) {
//     for (int i = 0; i < data->N; i++) {
//         // printf("i: %d (%d - %f)\n", i, data->r[i], f_hat->data[i]);
//         err = ((float) data->r[i]) - f_hat->data[i];
//         abs_error += fabs(err);
//         sq_error  += err * err;
//         if (fabs(err) < min_error) min_error = fabs(err);
//         if (fabs(err) > max_error) max_error = fabs(err);
//     }

//     printf("abs_error:    %lf\n", abs_error);
//     printf("sq_error:     %lf\n", sq_error);
//     printf("MSE:          %lf\n", sq_error / data->N);
//     printf("Average err:  %lf\n", abs_error / data->N);

//     ErrorReport report;

//     report.avg_error = abs_error / data->N;
//     report.max_error = max_error;
//     report.min_error = min_error;
//     report.mse = sq_error / data->N;
//     report.total_error = abs_error;

//     freePrecomp(precomp);
//     freeSphericalModel(model);
//     Matrix_free_f(f_hat);

//     return report;
// }