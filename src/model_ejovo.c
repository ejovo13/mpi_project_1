
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ejovo.h"
#include "ejovo/matrix/generic_macros.h"
#include "harmonics.h"
#include "gmp.h"

#define LMAX 50

/**========================================================================
 *!                           Looping macros
 *========================================================================**/
#define FOR_I_IN(mat, inner) for (size_t __i = 0, __n = Matrix_size_d(mat), i = 0; __i < __n; __i++, i = mat->data[__i]) { \
    i = mat->data[__i]; \
    inner \
    }
// iterate along the rows
#define FOR_I(mat, inner) for (size_t i = 0, __n = mat->nrows; i < __n; i++) { inner \ }
// #define FOR_N(mat, inner) for (size_t i = 0; i < Matrix_size_d(mat))
#define FOR_IJ(mat, inner) for (size_t i = 0; i < mat->nrows; i++) { \
    for (size_t j = 0; j < mat->ncols; j++) { \
        inner \
    }\
}

// iterate along the cols
#define FOR_J(mat, inner) for (size_t __j = 0, __n = mat->nrows)

iso_model *compute_model(int lmax, const char *filename, int n_theta, int n_phi);
iso_trained *train_model(iso_model *iso, int nrunsdouble, double alpha_zero);
iso_trained *train_model_both(iso_model *iso, int nruns, double alpha_zero);
void update_prediction(iso_model *iso);
void print_predictions(const iso_model *iso);
void report_training(const iso_trained *trained, const char *prefix);
void print_coeff(const iso_model *iso);
iso_trained *train_model_stochastic(iso_model *iso, int nruns, int __ndraws, double alpha_zero); 

int main(int argc, char *argv[]) {

    // Indices that will be helpful when converting between a linear index and a pair (l, m)


    // initialize rng
    ejovo_seed();

    // Start toying with the small data set
    // data_points_iso *data = load_data_points_iso("ETOPO1_small.csv", 64800);
    // int LMAX = 0;
    // if (argc >= 2) LMAX = 
    /**========================================================================
     *!                           DATA_FILES.csv
     *========================================================================**/
    const char filename_u[] = "ETOPO1_unit.csv";
    const char filename_t[] = "ETOPO1_tiny.csv";
    const char filename_s[] = "ETOPO1_small.csv";
    const char filename_m[] = "ETOPO1_med.csv";

    const int t_u = 6,     p_u = 12;
    const int t_t = 18,    p_t = 36;
    const int t_s = 180,   p_s = 360;
    const int t_m = 540,   p_m = 1080;

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
    // iso_model *iso = compute_model(lmax, filename_s, t_s, p_s);
    iso_model *iso = compute_model(lmax, filename_m, t_m, p_m);

    write_iso(iso->data, "iso.csv");

    printf("============== Predictions ===============\n");
    update_prediction(iso);
    // print_predictions(iso);
    // print_coeff(iso);

    printf("MSE(iso) : %lf\n", compute_mse(iso));
    printf("err(iso) : %lf\n", compute_average_error(iso)); 
    printf("ran with lmax := %d\n", lmax);

    /**========================================================================
     *!                           Gradient testing
     *========================================================================**/

    // int N_runs = 250;
    // double alpha = .01;
    // int N_draws = 50;

    // iso_model *iso_clone1 = copy_iso(iso);
    // iso_model *iso_clone2 = copy_iso(iso);

    // // iso_trained *trained = train_model(iso, N_runs, alpha);
    // // iso_trained *trained_stoch = train_model_stochastic(iso_clone1, N_runs, N_draws, alpha * 10);
    // iso_trained *trained = train_model(iso_clone2, N_runs, alpha);
    
    // report_training(trained_stoch, "test");
    // report_training(trained, "test");


    // printf("Initial error: \n");
    // print_el_d(Matrix_first_d(MSE));

    // printf("\nFinal error\n");
    // print_el_d(Matrix_last_d(MSE));
    // printf("\n");

    // double change_mse = Matrix_last_d(MSE) - Matrix_first_d(MSE);
    // printf("Total change in MSE after training: ");
    // if (change_mse > 0) {
    //     color_red();
    //     printf("%lf\n", change_mse);
    //     color_reset();
    // } else {
    //     printf("%lf\n", change_mse);
    // }

    // printf("rms: %lf\n", sqrt(Matrix_last_d(MSE)));

    // writeModel(iso->model, iso->data, "trained");

    // f_hat = modelPredictN(iso, iso->data->N);
    // reshape_d(f_hat, iso->data->p, iso->data->t);


    FILE *prediction_file = fopen("prediction.csv", "w");
    update_prediction(iso);

    for (int i = 0; i < iso->data->N; i++) {

        int i_th = data_i_th(iso->data, i);
        int i_ph = data_i_ph(iso->data, i);
        
        fprintf(prediction_file, "%lf\t%lf\t%lf\n", phiw_to_lamp(vecat_d(iso->data->ph, i_ph)),
                                                    thew_to_phip(vecat_d(iso->data->th, i_th)),
                                                    vecat_d(iso->f_hat, i));

    }

    fclose(prediction_file);

    // // Matrix_print_d(f_hat);

    // // Print the difference between f_hat and the observed data point
    // Matrix_d *diff = Matrix_subtract_d(f_hat, iso->data->r);

    // // Matrix_print_d(diff);

    // printf("Max diff: %lf\n", max_d(diff));
    // printf("Min diff: %lf\n", min_d(diff));
    // printf("Max_abs diff: %lf\n", maxabs_d(diff));
    // printf("std diff: %lf\n", std_d(diff));
    // printf("total_error: %lf\n", sumabs_d(diff));
    // printf("Alternative mse: %lf\n", mean_squared_d(diff));
    // printf("Meanabs: %lf\n", sumabs_d(diff) / Matrix_size_d(diff));
    // // printf("Average ? error: %lf\n", mean_d(diff));
    // // printf("model dt and dp: %lf, %lf\n", iso->data->dt, iso->data->dp);

    // // printf("mse(manual) : %lf\n", compute_mse(iso));


    // // compute_model(5, filename_s, 5, 10);
}


// Given an lmax, compute the values of C_lm and S_lm and return the values
// in a Vector of the form
//
//  Coeff = [ C_lm S_lm ]
//
// where the length of Coeff is 2 * (lmax + 1) * (lmax + 2) / 2
iso_model *compute_model(int lmax, const char *datafile, int n_theta, int n_phi) {

    const int ll = (lmax + 1) * (lmax + 2) / 2;

    printf("[compute_model] Constructing model { lmax: %d, t: %d, p: %d}\n", lmax, n_theta, n_phi);
    data_iso *data = load_data_iso(datafile, n_theta, n_phi);


    // use lmax to initialize data->l_indices 
    data->l_indices = Matrix_new_i(1, ll);
    data->m_indices = Matrix_new_i(1, ll);

    int i = 0;
    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            data->l_indices->data[i] = l;
            data->m_indices->data[i] = m;
            i++;
        }
    }

    iso_model *iso = newModel(data, lmax);
    writeModel(iso->model, data, "");    

    return iso;
}

void update_prediction(iso_model *iso) {
    Matrix_d *f_hat = modelPredictN(iso, iso->data->t * iso->data->p);
    // Matrix_d *f_hat = modelPredictN(iso, 1);
    reshape_d(f_hat, iso->data->p, iso->data->t);
    iso->f_hat = f_hat;
}

void print_predictions(const iso_model *iso) {


    for (int i = 0; i < iso->data->t * iso->data->p; i++) {

        int i_ph = i / iso->data->t;
        int i_th = i % iso->data->t;

        printf("Predictions stored in f_hat: %p\n", iso->f_hat);

        printf("th: %lf\t ph: %lf f(th, ph) = %lf\t\t = f(ph, lam) = ph: %lf\t lam: %lf\n", 
            vecat_d(iso->data->th, i_th), 
            vecat_d(iso->data->ph, i_ph), 
            vecat_d(iso->f_hat, i), 
            thew_to_phip(vecat_d(iso->data->th, i_th)),
            phiw_to_lamp(vecat_d(iso->data->ph, i_ph)));
    }
}

void print_coeff(const iso_model *iso) {
    printf("\n Clm \n");
    Matrix_print_d(iso->model->C_lm);
    printf("\n Slm \n");
    Matrix_print_d(iso->model->S_lm);
}

void report_training(const iso_trained *trained, const char *prefix) {

    const iso_model *iso = trained->iso;

    const char coeff[] = "coeff_series.csv";

    FILE *coeff_file = fopen(prepend(coeff, prefix), "w");

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
    for (int t = 0; t < trained->count; t++) {

        fprintf(coeff_file, "%d,", t);
        // iterate for all l,m in (0 .. L) X (0 .. l)
        for (size_t i = 0; i < iso->model->ll; i++) {
            fprintf(coeff_file, "%lf,", matat_d(trained->Clm_history, i, t));
        }

        for (size_t i = 0; i < iso->model->ll; i++) {
            fprintf(coeff_file, "%lf,", matat_d(trained->Slm_history, i, t));
        }
        fprintf(coeff_file, "%lf\n", trained->MSE->data[t]);
    }
}

// ndraws is the number of stochastic simulations to run. Make sure this number is greater than 0
iso_trained *train_model(iso_model *iso, int nruns, double alpha_zero) {

    printf("\nBeginning training...\n");

    // compute an initial MSE 
    double mse_initial = compute_mse(iso);

    const double eps = 0.00000005;

    double alpha = alpha_zero;
    // If our improvement affects less than this amount, stop training

    // so let's try and run the model 100 times to drop the MSE
    int count = 0;

    // Parameters during learning process
    Matrix_d *MSE = Matrix_new_d(1, nruns);
    Matrix_d *dMSE = Matrix_new_d(1, nruns - 1);
    Matrix_d *C_lm_t = Matrix_new_d(iso->model->ll, nruns);
    Matrix_d *S_lm_t = Matrix_new_d(iso->model->ll, nruns);
    Matrix_d *MSE_bench = Matrix_new_d(1, nruns);
    Matrix_d *grad_bench = Matrix_new_d(1, nruns);

    Clock *clock = Clock_new();

    for (int t = 0; t < nruns; t++) {

        // try and come up with a more robust learning rate parameter
        alpha = .99 * alpha;

        Clock_tic(clock);
            // Try computing the stochastic gradient
            Matrix_d *grad = compute_gradient(iso->data, iso->model);
            vecnormalize_d(grad);
        Clock_toc(clock);

        double grad_time = elapsed_time(clock);

        Clock_tic(clock);
            adjust_parameters(grad, iso->model, alpha);
        Clock_toc(clock);

        double adjustment_time = elapsed_time(clock);

        // Record the updated changes
        for (size_t j = 0; j < iso->model->ll; j++) {
            *matacc_d(C_lm_t, j, t) = vecat_d(iso->model->C_lm, j);
            *matacc_d(S_lm_t, j, t) = vecat_d(iso->model->S_lm, j);
        }

        Clock_tic(clock);
            double mse = compute_mse(iso);
        Clock_toc(clock);

        double mse_time = elapsed_time(clock);

        // Record timing variables
        MSE->data[t]        = mse;
        MSE_bench->data[t]  = mse_time;
        grad_bench->data[t] = grad_time; 

        count++;

        // report the change in mse:
        if (t > 0) {
            
            dMSE->data[t - 1] = MSE->data[t] - MSE->data[t - 1];

            // if (dMSE->data[t - 1] > 0) {
            //     // perror("Increasing MSE of mode\n");
            //     // exit(2);
            // }
            if (fabs(dMSE->data[t - 1]) < eps) {
                //
                // break;
                printf("Exiting for loop\n");
                break;
            }
        }


    }

    double mse_final = compute_mse(iso);
    

    // Matrix_d *MSE_econ = Matrix_new_d(1, count);
    // Matrix_d *dMSE_econ = Matrix_new_d(1, nruns - 1);
    Matrix_d *C_lm_t_econ = Matrix_new_d(iso->model->ll, count);
    Matrix_d *S_lm_t_econ = Matrix_new_d(iso->model->ll, count);
    Matrix_d *MSE_econ = take_d(MSE, count);
    Matrix_d *dMSE_econ = take_d(dMSE, count);

    assert(Matrix_size_d(MSE_econ) == (size_t) count);
    assert(Matrix_size_d(dMSE_econ) == (size_t) count);

    // Matrix_print_d(C_lm_t);
    // Matrix_print_d(S_lm_t);
            
    // Fill in the matrices only up until the count;
    FOR_IJ(C_lm_t_econ,
        *matacc_d(C_lm_t_econ, i, j) = matat_d(C_lm_t, i, j);
        *matacc_d(S_lm_t_econ, i, j) = matat_d(S_lm_t, i, j);
    )

    // allocate space for final trained model

    iso_trained *trained = (iso_trained *) malloc(sizeof(*trained));

    trained->count = count;
    trained->Clm_history = C_lm_t_econ;
    trained->Slm_history = S_lm_t_econ;
    trained->MSE = MSE_econ;
    trained->dMSE = dMSE_econ;
    trained->iso = iso;
    
    Matrix_free_d(C_lm_t);
    Matrix_free_d(S_lm_t);
    Matrix_free_d(MSE);
    Matrix_free_d(dMSE);

    printf("Training concluded, Matrices freed\n");
    printf("Initial MSE: %lf\n", mse_initial);
    printf("Final   MSE: %lf\n", mse_final);
    printf("trained %d iterations\n", trained->count);

    return trained;
}

// ndraws is the number of stochastic simulations to run. Make sure this number is greater than 0
iso_trained *train_model_stochastic(iso_model *iso, int nruns, int __ndraws, double alpha_zero) {

    const int ndraws = __ndraws > 0 ? __ndraws : 1;

    printf("\nBeginning training...\n");

    // compute an initial MSE 
    double mse_initial = compute_mse(iso);

    const double eps = 0.00000005;

    double alpha = alpha_zero;
    // If our improvement affects less than this amount, stop training

    // so let's try and run the model 100 times to drop the MSE
    int count = 0;

    // Parameters during learning process
    Matrix_d *MSE_hat = Matrix_new_d(1, nruns);
    Matrix_d *dMSE_hat = Matrix_new_d(1, nruns - 1);
    Matrix_d *C_lm_t = Matrix_new_d(iso->model->ll, nruns);
    Matrix_d *S_lm_t = Matrix_new_d(iso->model->ll, nruns);
    Matrix_d *MSE_hat_bench = Matrix_new_d(1, nruns);
    Matrix_d *grad_hat_bench = Matrix_new_d(1, nruns);

    Clock *clock = Clock_new();

    for (int t = 0; t < nruns; t++) {

        // try and come up with a more robust learning rate parameter
        alpha = .99 * alpha;

        Clock_tic(clock);
            // Try computing the stochastic gradient
            Matrix_d *grad_hat = compute_stochastic_gradient(iso, ndraws);
            vecnormalize_d(grad_hat);
        Clock_toc(clock);

        double grad_hat_time = elapsed_time(clock);

        Clock_tic(clock);
            adjust_parameters(grad_hat, iso->model, alpha);
        Clock_toc(clock);

        double adjustment_time = elapsed_time(clock);

        // Record the updated changes
        for (size_t j = 0; j < iso->model->ll; j++) {
            *matacc_d(C_lm_t, j, t) = vecat_d(iso->model->C_lm, j);
            *matacc_d(S_lm_t, j, t) = vecat_d(iso->model->S_lm, j);
        }

        Clock_tic(clock);
            double mse_hat = estimate_mse(iso, ndraws);
        Clock_toc(clock);

        double mse_hat_time = elapsed_time(clock);

        // Record timing variables
        MSE_hat->data[t]        = mse_hat;
        MSE_hat_bench->data[t]  = mse_hat_time;
        grad_hat_bench->data[t] = grad_hat_time; 

        count++;

        // report the change in mse:
        if (t > 0) {
            
            dMSE_hat->data[t - 1] = MSE_hat->data[t] - MSE_hat->data[t - 1];

            // if (dMSE_hat->data[t - 1] > 0) {
            //     perror("Increasing MSE of mode\n");
            //     // exit(2);
            // }
            if (fabs(dMSE_hat->data[t - 1]) < eps) {
                //
                // break;
                printf("Exiting for loop\n");
                break;
            }
        }


    }

    double mse_final = compute_mse(iso);
    

    // Matrix_d *MSE_econ = Matrix_new_d(1, count);
    // Matrix_d *dMSE_econ = Matrix_new_d(1, nruns - 1);
    Matrix_d *C_lm_t_econ = Matrix_new_d(iso->model->ll, count);
    Matrix_d *S_lm_t_econ = Matrix_new_d(iso->model->ll, count);
    Matrix_d *MSE_econ = take_d(MSE_hat, count);
    Matrix_d *dMSE_econ = take_d(dMSE_hat, count);

    assert(Matrix_size_d(MSE_econ) == (size_t) count);
    assert(Matrix_size_d(dMSE_econ) == (size_t) count);

    // Matrix_print_d(C_lm_t);
    // Matrix_print_d(S_lm_t);
            
    // Fill in the matrices only up until the count;
    FOR_IJ(C_lm_t_econ,
        *matacc_d(C_lm_t_econ, i, j) = matat_d(C_lm_t, i, j);
        *matacc_d(S_lm_t_econ, i, j) = matat_d(S_lm_t, i, j);
    )

    // allocate space for final trained model

    iso_trained *trained = (iso_trained *) malloc(sizeof(*trained));

    trained->count = count;
    trained->Clm_history = C_lm_t_econ;
    trained->Slm_history = S_lm_t_econ;
    trained->MSE = MSE_econ;
    trained->dMSE = dMSE_econ;
    trained->iso = iso;
    
    Matrix_free_d(C_lm_t);
    Matrix_free_d(S_lm_t);
    Matrix_free_d(MSE_hat);
    Matrix_free_d(dMSE_hat);

    printf("Training concluded, Matrices freed\n");
    printf("Initial MSE: %lf\n", mse_initial);
    printf("Final   MSE: %lf\n", mse_final);
    printf("trained %d iterations\n", trained->count);

    return trained;
}

iso_trained *train_model_both(iso_model *iso, int nruns, double alpha_zero) {

    printf("\nBeginning training...\n");

    const double eps = 0.00000005;

    double alpha = alpha_zero;
    // If our improvement affects less than this amount, stop training

    // so let's try and run the model 100 times to drop the MSE
    int count = 0;

    // Parameters during learning process
    Matrix_d *grad = NULL;
    Matrix_d *MSE = Matrix_new_d(1, nruns);
    Matrix_d *dMSE = Matrix_new_d(1, nruns - 1);
    Matrix_d *C_lm_t = Matrix_new_d(iso->model->ll, nruns);
    Matrix_d *S_lm_t = Matrix_new_d(iso->model->ll, nruns);

    // Temporarily find out how much different a stochastic approximation is
    Matrix_d *grad_distance = Matrix_new_d(1, nruns);

    Clock *clock = Clock_new();

    for (int t = 0; t < nruns; t++) {

        // I would also like to store the time series of Clm [Slm] for every iteration
        // to show how learning has changed the MSE as well as the actual coefficients 
        // themselves

        alpha = .99 * alpha;

        printf("Computing gradient with alpha = %lf and t = %d\n", alpha, t);

        Clock_tic(clock);
        // Try computing the stochastic gradient
        Matrix_d *stochastic_grad = compute_stochastic_gradient(iso, iso->data->N / 5);
        grad = compute_gradient(iso->data, iso->model);

        vecnormalize_d(stochastic_grad);
        vecnormalize_d(grad);

        // Difference between the two gradients
        Matrix_d *grad_diff = Matrix_subtract_d(grad, stochastic_grad);
        double grad_dist = vecnorm_d(grad_diff);
        printf("|| grad_diff || = %lf\n", grad_dist);       

        grad_distance->data[t] = grad_dist;

        // Compare the gradient and the stochastic gradient

        printf("Gradient vector: ");
        Matrix_print_d(grad);
        printf("Stochastic Gradient: ");
        Matrix_print_d(stochastic_grad);

        Clock_toc(clock);
        adjust_parameters(grad, iso->model, alpha);

        printf("normalized gradient in %lfs...\n", elapsed_time(clock));
        // have alpha decrease

        // Record the updated changes
        for (size_t j = 0; j < iso->model->ll; j++) {
            *matacc_d(C_lm_t, j, t) = vecat_d(iso->model->C_lm, j);
            *matacc_d(S_lm_t, j, t) = vecat_d(iso->model->S_lm, j);
        }


        Clock_tic(clock);
        MSE->data[t] = compute_mse(iso);
        Clock_toc(clock);

        double mse_time = elapsed_time(clock);

        Clock_tic(clock);
        double mse_hat = estimate_mse(iso, iso->data->N / 10);
        // double mse_hat = estimate_mse(iso, iso->data->N);
        Clock_toc(clock);

        double mse_hat_time = elapsed_time(clock);


        printf("Computed mse in %lfs...\n", mse_time);
        printf("Estimated mse: %lf with mse_hat: %lf in %lfs, mse - mse_hat = %lf\n", MSE->data[t], mse_hat, mse_hat_time, mse_hat - MSE->data[t]);
        count++;

        // report the change in mse:
        if (t > 0) {
            
            dMSE->data[t - 1] = MSE->data[t] - MSE->data[t - 1];

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

    

    // Matrix_d *MSE_econ = Matrix_new_d(1, count);
    // Matrix_d *dMSE_econ = Matrix_new_d(1, nruns - 1);
    Matrix_d *C_lm_t_econ = Matrix_new_d(iso->model->ll, count);
    Matrix_d *S_lm_t_econ = Matrix_new_d(iso->model->ll, count);
    Matrix_d *MSE_econ = take_d(MSE, count);
    Matrix_d *dMSE_econ = take_d(dMSE, count);

    assert(Matrix_size_d(MSE_econ) == (size_t) count);
    assert(Matrix_size_d(dMSE_econ) == (size_t) count);

    // Matrix_print_d(C_lm_t);
    // Matrix_print_d(S_lm_t);
            
    // Fill in the matrices only up until the count;
    FOR_IJ(C_lm_t_econ,
        *matacc_d(C_lm_t_econ, i, j) = matat_d(C_lm_t, i, j);
        *matacc_d(S_lm_t_econ, i, j) = matat_d(S_lm_t, i, j);
    )

    // allocate space for final trained model

    iso_trained *trained = (iso_trained *) malloc(sizeof(*trained));

    trained->count = count;
    trained->Clm_history = C_lm_t_econ;
    trained->Slm_history = S_lm_t_econ;
    trained->MSE = MSE_econ;
    trained->dMSE = dMSE_econ;
    trained->iso = iso;
    
    Matrix_free_d(C_lm_t);
    Matrix_free_d(S_lm_t);

    return trained;
}
