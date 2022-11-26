#include "geodesy.h"

/**========================================================================
 * ?                          learn.c
 * @brief   : Functions dealing with learning and predicting
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-16
 *========================================================================**/


// Predict any single arbitrary point f(\theta, \phi)
// Runs in O(l^2) time 
double modelPredict(const iso_model *iso, double theta, double phi) {

    // Compute the P values
    Matrix_d *P_lm = Matrix_new_d(1, iso->model->ll);

    // computeP(iso->coeff, P_lm->data, cos(theta));
    computeP(iso->coeff, P_lm->data, cos(theta));

    // Now that we have plm, compute the sum using the model's coefficients
    double sum = 0;

    for (size_t l = 0; l <= iso->model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {
            // Plm(cos\theta)*[clm cos(m\phi) + slm sin(m\phi)]
            sum += vecat_d(P_lm, PT(l, m)) * (model_C(iso->model, l, m) * cos(m * phi) +
                                              model_S(iso->model, l, m) * sin(m * phi));

            // printf("Summing Clm: %lf and Slm: %lf\n", model_C(iso->model, l, m), model_S(iso->model, l, m));
        }
    }

    return sum;
}

// Predict a single point that belongs to the data set
double modelPredictDataPoint(const iso_model *iso, int i) {

    const Matrix_d *P_lm_th = iso->model->P_lm_th;

    int i_ph = data_i_ph(iso->data, i);
    int i_th = data_i_th(iso->data, i);

    double ph = vecat_d(iso->data->th, i_ph);
    // double th = vecat_d(iso->data->ph, i_th);

    double sum = 0;

    for (size_t l = 0; l <= iso->model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {
            // Plm(cos\theta)*[clm cos(m\phi) + slm sin(m\phi)]
            sum += matat_d(P_lm_th, i_th, PT(l, m)) * (model_C(iso->model, l, m) * cos(m * ph)) +
                                                       model_S(iso->model, l, m) * sin(m * ph);
        }
    }

    return sum;
}

Matrix_d *modelPredictDataPoints(const iso_model *iso, Matrix_i *indices) {

    const int n = Matrix_size_i(indices);
    Matrix_d *predictions = Matrix_new_d(1, n);

    for (int i = 0; i < n; i++) {
        predictions->data[i] = modelPredictDataPoint(iso, indices->data[i]);
    }

    return predictions;
}

// Convert a data set to an output predicted by the model
Matrix_d *modelPredictData(const iso_model *iso) {

    Matrix_d *m = Matrix_new_d(iso->data->p, iso->data->t);

    // For all of the data points
    for (int i = 0; i < iso->data->N; i++) {
        m->data[i] = modelPredictDataPoint(iso, i);
    }

    return m;
}

// Compute the first n predictions of models data
Matrix_d *modelPredictN(const iso_model *iso, int __n) {

    int n = iso->data->N < __n ? iso->data->N : __n;

    const Matrix_d *th = iso->data->th;
    const Matrix_d *ph = iso->data->ph;

    Matrix_d *predictions = Matrix_new_d(1, n);

    for (int i = 0; i < n; i++) {

        int i_ph = i / iso->data->t;
        int i_th = i % iso->data->t;

        predictions->data[i] = modelPredict(iso, vecat_d(th, i_th), vecat_d(ph, i_ph));
        // printf("=======================\n");
    }

    return predictions;
}

// Compute gradient in the form 
// [dMSE/dc00 dMSE/dc10 ... dMSE/dclm
//  dMSE/ds00 dMSE/ds10 ... dMSE/dslm]
// Assume model has C_lm and S_lm and P_lm_th
Matrix_d *compute_gradient(const data_iso *data, const spherical_model *model) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    const Matrix_d *clm = NULL;
    const Matrix_d *slm = NULL;
    const Matrix_d *pcs = NULL;

    // printf("Hi im in gradient\n");

    // printf("%p, %p, %p\n", model->clm, model->slm, model->pcs);

    // check if the model has the models already computed
    // if (model->clm != NULL && model->slm != NULL && model->pcs != NULL) {
    if (model->clm != NULL && model->slm != NULL) {
    // We will actually need to recompute the indices for a general purpose calculation
    // if (false && model->clm != NULL && model->slm != NULL && model->pcs != NULL) {
        // printf("Already computed coefficients\n");
        clm = model->clm; // 2 x ll matrix
        slm  = model->slm;
        // pcs = model->pcs; // sum(clm * Plmcos * cos)_i for i in 1..N
    } else {
        clm = compute_mse_coeff_clm(data, model); // 1 x ll matrix
        slm  = compute_mse_coeff_slm(data, model); // 1 x ll matrix
    }

    // Unless I am mistaken, we HAVE to recompute pcs
    pcs = compute_pcs(data, model, clm, slm); // sum(clm * Plmcos * cos)_i for i in 1..N
    // allocate space for output array

    // double *grad = (double *) malloc(sizeof(*grad) * 2 * ll);
    Matrix_d *grad = Matrix_new_d(2, ll);

    for (size_t l = 0; l <= model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {

            grad->data[PT(l, m)] = compute_gradient_clm(data, clm, pcs, l, m);
            grad->data[PT(l, m) + ll] = compute_gradient_slm(data, slm, pcs, l, m);

        }
    }

    return grad;
}

Matrix_d *compute_stochastic_gradient(const iso_model *iso, int n) {


    const int ll = iso->model->ll;

    Matrix_i *indices = runif_i(n, 0, iso->data->N);
    // printf("Using random indices: \n");
    // Matrix_print_i(indices);

    const Matrix_d *clm = iso->model->clm; // 2 x ll matrix
    const Matrix_d *slm = iso->model->slm;

    if (iso->model->clm != NULL && iso->model->slm != NULL) {
    // We will actually need to recompute the indices for a general purpose calculation
        clm = iso->model->clm; // 2 x ll matrix
        slm  = iso->model->slm;
    } else {
        clm = compute_mse_coeff_clm(iso->data, iso->model); // 1 x ll matrix
        slm  = compute_mse_coeff_slm(iso->data, iso->model); // 1 x ll matrix
    }

    const Matrix_d *pcs = compute_pcs_indices(iso, indices);

    // generate random indices

    // printf("Hi im in gradient\n");

    // printf("%p, %p, %p\n", iso->model->clm, iso->model->slm, iso->model->pcs);
    Matrix_d *grad = Matrix_new_d(2, ll);

    for (size_t l = 0; l <= iso->model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {

            grad->data[PT(l, m)]      = compute_gradient_clm_points(iso, clm, pcs, l, m, indices);
            grad->data[PT(l, m) + ll] = compute_gradient_slm_points(iso, slm, pcs, l, m, indices);

        }
    }

    Matrix_free_i(indices);

    return grad;

}

// Assume that a model object has been fully initialized
double compute_mse(const iso_model *iso) {

    // const spherical_model *model = iso->model;
    const data_iso        *data  = iso->data;


    // const int ll = model->ll;

    double sum = 0;

    // Compute MSE by getting a list of predictions
    // const Matrix_d *predictions = modelPredictN(iso, iso->data->N);
    // const Matrix_d *predictions = modelPredictData(iso);
    const Matrix_d *predictions = modelPredictData(iso);

    // now compute the difference between the two
    for (int i = 0; i < data->N; i++) {
        double diff = vecat_d(predictions, i) - data->r[i];
        sum += diff * diff;
    }

    return sum / data->N;
}

double estimate_mse(const iso_model *iso, int n) {

    // estimate the MSE using only i data points
    Matrix_i *indices = runif_i(n, 0, iso->data->N);

    // printf("Sampling from indices: \n");
    // Matrix_print_i(indices);

    // const spherical_model *model = iso->model;
    const data_iso        *data  = iso->data;

    double sum = 0;

    const Matrix_d *predictions = modelPredictDataPoints(iso, indices);

    for (int i = 0; i < n; i++) {
        int index = indices->data[i] - 1;
        double diff = vecat_d(predictions, i) - data->r[index];
        sum += diff * diff;
    }

    return sum / n;
}

double compute_average_error(const iso_model *iso) {

    // const spherical_model *model = iso->model;
    const data_iso        *data  = iso->data;


    // const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    double sum = 0;

    // Compute MSE by getting a list of predictions
    const Matrix_d *predictions = modelPredictN(iso, iso->data->N);

    // now compute the difference between the two
    for (int i = 0; i < data->N; i++) {
        sum += vecat_d(predictions, i) - data->r[i];
    }

    return sum / data->N;
}

// Given a gradient matrix and a model, tweak the values of 
// C_lm to improve the score
// gradient has size (N * 2) x ((lmax + 1)(lmax + 2)/2)
void adjust_parameters(const Matrix_d *grad, spherical_model *model, double alpha) {

    // const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    // create a scaled version of the gradient vector (not necessary, but since the gradient vector should 
    // be reasonably small, we can perform this operation with very little cost)
    Matrix_d *grad_scaled = Matrix_clone_d(grad); // FAST memory copy
    matmultscalar_d(grad_scaled, -alpha); // Just multiply the first row

    MatIter_row_add_row_d(Matrix_row_begin_d(model->C_lm, 0), 
                          Matrix_row_end_d(model->C_lm, 0), 
                          Matrix_row_begin_d(grad_scaled, 0));

    MatIter_row_add_row_d(Matrix_row_begin_d(model->S_lm, 0),
                          Matrix_row_end_d(model->S_lm, 0),
                          Matrix_row_begin_d(grad_scaled, 1));
   
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

        // double adjustment_time = elapsed_time(clock);

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

        // double adjustment_time = elapsed_time(clock);

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