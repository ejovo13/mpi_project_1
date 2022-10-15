
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ejovo.h"
#include "harmonics.h"
#include "gmp.h"

#define LMAX 50

iso_model *compute_model(int lmax, const char *filename, int n_theta, int n_phi);

// const int nrows = (LMAX + 1) * (LMAX + 2) / 2;

int main(int argc, char *argv[]) {


    // Start toying with the small data set
    // data_points_iso *data = load_data_points_iso("ETOPO1_small.csv", 64800);
    // int LMAX = 0;
    // if (argc >= 2) LMAX = 
    /**========================================================================
     *!                           ETOPO1_unit.csv
     *========================================================================**/
    const char filename_u[] = "ETOPO1_unit.csv";
    const int t_u = 6;
    const int p_u = 12;

    /**========================================================================
     *!                           ETOPO1_tiny.csv
     *========================================================================**/
    const char filename_t[] = "ETOPO1_tiny.csv";
    const int t_t = 18;
    const int p_t = 36;

    /**========================================================================
     *!                           ETOPO1_small.csv
     *========================================================================**/
    const char filename_s[] = "ETOPO1_small.csv";
    const int t_s = 180;
    const int p_s = 360;

    /**========================================================================
     *!                           ETOPO1_med.csv
     *========================================================================**/
    const char filename_m[] = "ETOPO1_med.csv";
    const int t_m = 540;
    const int p_m = 1080;

    // get lmax from the argument pasaed
    int lmax = 0;

    if (argc >= 2) {
        sscanf(argv[1], "%d", &lmax); 
        printf("lmax set to %d with %d args\n", lmax, argc);
    } else {
        fprintf(stderr, "Please enter a value for lmax\n ");
        fprintf(stderr, "ex:  %s 10\n", argv[0]);
        exit(2);
    }

    // iso_model *iso = compute_model(5, filename_s, 4, 5);
    // iso_model *iso = compute_model(lmax, filename_u, t_u, p_u);
    // iso_model *iso = compute_model(lmax, "ETOPO1_unif.csv", t_u, p_u);
    // iso_model *iso = compute_model(lmax, filename_t, t_t, p_t);
    // iso_model *iso = compute_model(lmax, filename_t, t_t, p_t);

    iso_model *iso = compute_model(lmax, filename_s, t_s, p_s);

    write_iso(iso->data, "iso.csv");

    // int n = 10;

    // print phi
    // Matrix_print_d(iso->data->ph);
    // Matrix_print_d(map_d(iso->data->ph, phip_to_thew));

    Matrix_d *copy = Matrix_clone_d(iso->data->ph);
    apply_d(copy, phiw_to_lamp);
    Matrix_print_d(copy);

    printf("============== Predictions ===============\n");
    Matrix_d *f_hat = modelPredictN(iso, t_u * p_u);
    // Matrix_d *f_hat = modelPredictN(iso, 1);
    reshape_d(f_hat, p_u, t_u);
    for (int i = 0; i < t_u * p_u; i++) {

        int i_ph = i / iso->data->t;
        int i_th = i % iso->data->t;

        printf("th: %lf\t ph: %lf f(th, ph) = %lf\t\t = f(ph, lam) = ph: %lf\t lam: %lf\n", 
            vecat_d(iso->data->th, i_th), 
            vecat_d(iso->data->ph, i_ph), 
            vecat_d(f_hat, i), 
            thew_to_phip(vecat_d(iso->data->th, i_th)),
            phiw_to_lamp(vecat_d(iso->data->ph, i_ph)));
            // vecat_d(iso->data->ph, i_ph));
            // 5.0);
            // vecat_d(iso->data->ph, i_ph));

        // printf("%lf\t%lf\t%lf\n", vecat_d(iso->data->th, i_th) - )
    }

    // So I should have access to the full list of predicitons and the full list 
    // of f
    // Matrix_print_d(iso->data->r);

    // Matrix_print_d(iso->model->C_lm);
    // Matrix_print_d(iso->model->S_lm);

    printf("MSE(iso) : %lf\n", compute_mse(iso));
    printf("err(iso) : %lf\n", compute_average_error(iso));

    // printf("\t========================f_hat ====================\t\n"); Matrix_print_d(f_hat);

    // Need a function to output the prediction of the data to a file in NASA's form


    
    // Coefficients of Clm and Slm
    // printf("\t======================= C_lm := ===================\t\n");
    // Matrix_print_d(iso->model->C_lm);
    // printf("\t======================= S_lm := ===================\t\n");
    // Matrix_print_d(iso->model->S_lm);

    // printf("P_lm_th: ");
    // Matrix_print_d(iso->model->P_lm_th);

    printf("ran with lmax := %d\n", lmax);
    // printf("rms(ise) : %lf\n", sqrt(compute_mse(iso)));

    // Manually set the model's coefficients
    // Matrix_fill_d(iso->model->C_lm, 0);
    // Matrix_fill_d(iso->model->S_lm, 0);
    // iso->model->C_lm->data[0] = 1000; // set the first coefficient to 1000, then generate a model


    // f_hat = modelPredictN(iso, t_u * p_u);
    // Matrix_d *f_hat = modelPredictN(iso, 1);
    // reshape_d(f_hat, p_u, t_u);
    // Matrix_print_d(f_hat);



    /**========================================================================
     *!                           Gradient testing
     *========================================================================**/
    int N_runs = 5;
    double alpha = 100;

    printf("\nBeginning training...\n");

    const double eps = 0.00000005;
    // If our improvement affects less than this amount, stop training

    // so let's try and run the model 100 times to drop the MSE
    int count = 0;

    // Parameters during learning process
    Matrix_d *grad = NULL;
    Matrix_d *MSE = Matrix_new_d(1, N_runs);
    Matrix_d *dMSE = Matrix_new_d(1, N_runs - 1);
    Matrix_d *C_lm_t = Matrix_new_d(iso->model->ll, N_runs);
    Matrix_d *S_lm_t = Matrix_new_d(iso->model->ll, N_runs);

    Clock *clock = Clock_new();

    for (int t = 0; t < N_runs; t++) {

        // I would also like to store the time series of Clm [Slm] for every iteration
        // to show how learning has changed the MSE as well as the actual coefficients 
        // themselves

        alpha = .99 * alpha;

        printf("Computing gradient with alpha = %lf and t = %d\n", alpha, t);

        Clock_tic(clock);
        // Try computing the stochastic gradient
        // grad = compute_stochastic_gradient(iso, 1000);
        grad = compute_gradient(iso->data, iso->model);
        vecnormalize_d(grad);
        Clock_toc(clock);
        adjust_parameters(grad, iso->model, alpha);

        printf("normalized gradient in %lfs...\n", elapsed_time(clock));
        // have alpha decrease



        Clock_tic(clock);
        MSE->data[t] = compute_mse(iso);
        Clock_toc(clock);

        double mse_time = elapsed_time(clock);

        Clock_tic(clock);
        // double mse_hat = estimate_mse(iso, 10000);
        double mse_hat = estimate_mse(iso, iso->data->N);
        Clock_toc(clock);

        double mse_hat_time = elapsed_time(clock);


        printf("Computed mse in %lfs...\n", mse_time);
        printf("Estimated mse: %lf with mse_hat: %lf in %lfs, mse - mse_hat = %lf\n", MSE->data[t], mse_hat, mse_hat_time, mse_hat - MSE->data[t]);

        // report the change in mse:
        if (t > 0) {
            
            dMSE->data[t - 1] = MSE->data[t] - MSE->data[t - 1];
            count ++;

            if (dMSE->data[t - 1] > 0) {
                perror("Increasing MSE of mode\n");
                // exit(2);
            }
            if (fabs(dMSE->data[t - 1]) < eps) {
                //
                // break;
                printf("Exiting for loop\n");
                break;
            }
        }
    }

    // For 1 in count, write some of the observed data to a file

    const char coeff[] = "coeff_series.csv";

    FILE *coeff_file = fopen(coeff, "w");

    // write header
    // for (size_t i = 0; i < iso->model->ll; i++) {    
    fprintf(coeff_file, "t,");

    for (size_t l = 0; l <= iso->model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {

            // write clm or slm
            fprintf(coeff_file, "c%lu_%lu,", l, m);
        }
    }

    for (size_t l = 0; l <= iso->model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {

            // write clm or slm
            fprintf(coeff_file, "s%lu_%lu,", l, m);
        }
    }
    fprintf(coeff_file, "mse\n");
    // fprintf(coeff_file, "\n");
    
    // Now write the actual data
    for (int t = 0; t < count; t++) {


        fprintf(coeff_file, "%d,", t);
        // iterate for all l,m in (0 .. L) X (0 .. l)
        for (size_t i = 0; i < iso->model->ll; i++) {
            fprintf(coeff_file, "%lf,", matat_d(C_lm_t, i, t));
        }

        for (size_t i = 0; i < iso->model->ll; i++) {
            fprintf(coeff_file, "%lf,", matat_d(C_lm_t, i, t));
        }
        fprintf(coeff_file, "%lf\n", MSE->data[t]);
    }

    

    // printf("o")
    reshape_d(dMSE, 1, count);
    reshape_d(MSE, 1, count);

    
    // // Matrix_d *diff
    // Matrix_print_d(dMSE);

    // printf("count :%d\n", count);

    // printf("\nCoefficients after gradient descent\n");
    // Matrix_print_d(iso->model->C_lm);
    // Matrix_print_d(iso->model->S_lm);
    // printf("\n");

    printf("\nMSE:\n");
    // Print the status of MSE
    Matrix_print_d(MSE);

    printf("Initial error: \n");
    print_el_d(Matrix_first_d(MSE));

    printf("\nFinal error\n");
    print_el_d(Matrix_last_d(MSE));
    printf("\n");

    double change_mse = Matrix_last_d(MSE) - Matrix_first_d(MSE);
    printf("Total change in MSE after training: ");
    if (change_mse > 0) {
        color_red();
        printf("%lf\n", change_mse);
        color_reset();
    } else {
        printf("%lf\n", change_mse);
    }

    printf("rms: %lf\n", sqrt(Matrix_last_d(MSE)));

    writeModel(iso->model, iso->data, "trained");

    f_hat = modelPredictN(iso, iso->data->N);
    reshape_d(f_hat, iso->data->p, iso->data->t);


    FILE *prediction_file = fopen("prediction.csv", "w");

    for (int i = 0; i < iso->data->N; i++) {

        int i_th = data_i_th(iso->data, i);
        int i_ph = data_i_ph(iso->data, i);
        
        fprintf(prediction_file, "%lf\t%lf\t%lf\n", phiw_to_lamp(vecat_d(iso->data->ph, i_ph)),
                                                    thew_to_phip(vecat_d(iso->data->th, i_th)),
                                                    vecat_d(f_hat, i));

    }

    fclose(prediction_file);

    // Matrix_print_d(f_hat);

    // Print the difference between f_hat and the observed data point
    Matrix_d *diff = Matrix_subtract_d(f_hat, iso->data->r);

    // Matrix_print_d(diff);

    printf("Max diff: %lf\n", max_d(diff));
    printf("Min diff: %lf\n", min_d(diff));
    printf("Max_abs diff: %lf\n", maxabs_d(diff));
    printf("std diff: %lf\n", std_d(diff));
    printf("total_error: %lf\n", sumabs_d(diff));
    printf("Alternative mse: %lf\n", mean_squared_d(diff));
    printf("Meanabs: %lf\n", sumabs_d(diff) / Matrix_size_d(diff));
    // printf("Average ? error: %lf\n", mean_d(diff));
    // printf("model dt and dp: %lf, %lf\n", iso->data->dt, iso->data->dp);

    // printf("mse(manual) : %lf\n", compute_mse(iso));


    // compute_model(5, filename_s, 5, 10);
}


// Given an lmax, compute the values of C_lm and S_lm and return the values
// in a Vector of the form
//
//  Coeff = [ C_lm S_lm ]
//
// where the length of Coeff is 2 * (lmax + 1) * (lmax + 2) / 2
iso_model *compute_model(int lmax, const char *filename, int n_theta, int n_phi) {

    printf("[compute_model] Constructing model { lmax: %d, t: %d, p: %d}\n", lmax, n_theta, n_phi);
    data_iso *data = load_data_iso(filename, n_theta, n_phi);
    iso_model *iso = newModel(data, lmax);
    writeModel(iso->model, data, "");    

    return iso;
}
