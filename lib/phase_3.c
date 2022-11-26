#include "geodesy.h"

// A prediction uses floats (to reduce size by a factor of 2) and is simple a length N row 
// vector that is computed from a model, precomputed values, and a data set. A prediction 
// will be able to be stored efficiently in binary so that we can effectively compute the 
// MSE associated 
Matrix_f *compute_prediction(const SphericalModel *model, const Precomp *precomp, const data_iso* data) {

    // A prediction is simply the evaluation of a Laplace Series evaluated at all the data points
    const int N = data->N;

    Matrix_f *f_hat = Matrix_new_f(1, N);

    for (int i = 0; i < N; i++) {
        *vecacc_f(f_hat, i) = compute_prediction_point(model, precomp, data, i);
    }

    return f_hat;
}

Matrix_f *compute_prediction_omp(const SphericalModel *model, const Precomp *precomp, const data_iso* data) {

        // A prediction is simply the evaluation of a Laplace Series evaluated at all the data points
    const int N = data->N;

    Matrix_f *f_hat = Matrix_new_f(1, N);

    #pragma omp parallel shared(f_hat) 
    {
        #pragma omp for
        for (int i = 0; i < N; i++) {
            *vecacc_f(f_hat, i) = compute_prediction_point(model, precomp, data, i);
        }

    }

    return f_hat;

}

Matrix_f *compute_prediction_mpi(const SphericalModel *model, const Precomp *precomp, const data_iso* data, const char * dataset_size, int world_size, int this_rank) {

    const int N = data->N;

    Matrix_f *f_hat = NULL;

    if (this_rank == 0) {
        f_hat = Matrix_new_f(1, N);
    } else {
        f_hat = Matrix_new_f(0, 0);
    }

    // printf("[%d] f_hat is a matrix of size: %d\n", this_rank, Matrix_size_f(f_hat));

    int this_workload = compute_workload(N, world_size, this_rank);
    Matrix_f *this_f_hat = Matrix_new_f(1, this_workload);

    assert(Matrix_size_f(this_f_hat) == this_workload);

    Matrix_i *start_end = compute_startend_array(N, world_size);

    // MPI_ONCE(
    //     printf("start end matrix: \n");
    //     Matrix_print_i(start_end);
    // )
    Matrix_i *displacements = compute_displacements(N, world_size);
    Matrix_i *recvcount = compute_workload_array(N, world_size);

    // MPI_ONCE(
    //     Matrix_print_i(displacements);
    //     Matrix_print_i(recvcount);
    // )
    // We want 
    int start = start_end->data[this_rank] + 1;
    int end   = start_end->data[this_rank + 1];

    // MPI_ORDERED(
    //     printf("[%d] predicting i in [%d, %d)\n", this_rank, start, end);
    //     printf("Model address: %p\n", model);
    // )

    // Vector_print_head_d(model->C_lm, 10);

    // MPI_ORDERED(
    //     printf("[%d] length of this_f_hat: %lu\n", this_rank, Matrix_size_f(this_f_hat));
    // )

    assert(data);
    assert(data->th);

    for (int i = start, count = 0; i <= end; i++, count++) {
        *vecacc_f(this_f_hat, count) = compute_prediction_point(model, precomp, data, i);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // MPI_ORDERED(
    //     printf("this_f_hat filled\n");
    // )
    
    MPI_Gatherv(
        this_f_hat->data,
        this_workload,
        MPI_FLOAT,
        f_hat->data,
        recvcount->data,
        displacements->data,
        MPI_FLOAT,
        0,
        MPI_COMM_WORLD
    );

    // Now Load the predictions to sent it to everyone
    
    MPI_Barrier(MPI_COMM_WORLD);

    

    char binary_in[100] = {0};
    sprintf(binary_in, "fhat_%s_%lu.bin", dataset_size, model->lmax);

    if (this_rank == 0) {
        printf("[compute_prediction_mpi] Saving predictions in %s\n", binary_in);
        save_prediction(f_hat, binary_in);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (this_rank != 0) 
        f_hat = load_prediction(binary_in, N);

    return f_hat; 

}

float compute_prediction_point(const SphericalModel *model, const Precomp *precomp, const data_iso *data, int i) {

    const Matrix_d *P_lm_th = precomp->Plm_th;

    int i_ph = data_i_ph(data, i);
    int i_th = data_i_th(data, i);

    double sum = 0;



    for (size_t l = 0; l <= model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {
            // Plm(cos\theta)*[clm cos(m\phi) + slm sin(m\phi)]
            sum += matat_d(P_lm_th, i_th, PT(l, m)) * (vecat_d(model->C_lm, PT(l, m)) * matat_d(precomp->cosmph, m, i_ph) +
                                                       vecat_d(model->S_lm, PT(l, m)) * matat_d(precomp->sinmph, m, i_ph));
        }
    }

    // printf("predicted value: %f\n", sum);
    return sum;
}

float compute_prediction_point_mpi(const SphericalModel *model, const Precomp *precomp, const data_iso *data, int i) {

    const Matrix_d *P_lm_th = precomp->Plm_th;

    int i_ph = data_i_ph(data, i);
    int i_th = data_i_th(data, i);

    double sum = 0;

    for (size_t l = 0; l <= model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {
            // Plm(cos\theta)*[clm cos(m\phi) + slm sin(m\phi)]
            sum += matat_d(P_lm_th, i_th, PT(l, m)) * (vecat_d(model->C_lm, PT(l, m)) * matat_d(precomp->cosmph, m, i_ph) +
                                                       vecat_d(model->S_lm, PT(l, m)) * matat_d(precomp->sinmph, m, i_ph));
        }
    }

    // printf("predicted value: %f\n", sum);
    return sum;
}



void save_prediction(const Matrix_f *f_hat, const char *binary_out) {

    FILE *bin = fopen(binary_out, "wb");
    if (bin == NULL) errx(1, "Binary prediction file does not exist");

    // Otherwise, write all N bytes to the binary file
    const size_t n_bytes = sizeof(float) * Matrix_size_f(f_hat);

    fwrite(f_hat->data, n_bytes, 1, bin);

    fclose(bin);

}

Matrix_f *load_prediction(const char *binary_in, int N) {

    if (N <= 0) errx(1, "N [%d] must be greater than 0\n", N);
    Matrix_f *f_hat = Matrix_new_f(1, N);

    FILE *bin = fopen(binary_in, "rb");
    if (bin == NULL) errx(1, "Binary prediction file does not exist");    

    const size_t n_bytes = sizeof(float) * N;


    fread(f_hat->data, n_bytes, 1, bin);

    fclose(bin);

    return f_hat;

}

// Driver segment of code to load predictions and write differences to files
void predict_stuff(const args_t *args, const SphericalModel *model, const Precomp *precomp) {
        // compute the predicted values and store them in a binary file
    // fhat_<size_dataset>_<lmodel>.bin
    if (args->predict) {

        char binary_prediction_out[100] = {0};
        sprintf(binary_prediction_out, "fhat_%s_%d.bin", args->size_dataset, args->lmodel);
        Matrix_f *f_hat = NULL;
        Clock *clock = Clock_new();

        FILE *test_open = fopen(binary_prediction_out, "rb");
        if (test_open == NULL || args->recompute) {

            /**========================================================================
             *!                           Compute Prediction
            *========================================================================**/
            printf("[main] Computing altitude predictions...\n");
            Clock_tic(clock);
            f_hat = compute_prediction(model, precomp, args->data);
            Clock_toc(clock);
            printf("[main] Time to compute prediction: %lfs\n", elapsed_time(clock));
            save_prediction(f_hat, binary_prediction_out);
            Vector_print_head_f(f_hat, 10);

            Matrix_free_f(f_hat);

        } else {

            printf("[main] Predictions already computed\n");
            fclose(test_open);
        }

        head_data(args->data);

        // We can't reach this point before having the predictions already computed
        if (args->diff) {

            char diff_out[100] = {0};
            sprintf(diff_out, "diff_%s_%d.csv", args->size_dataset, args->lmodel);

            // test_open = fopen(diff_out, "r");
            // if (test_open == NULL) {

                Clock_tic(clock);
                f_hat = load_prediction(binary_prediction_out, args->data->N);
                Clock_toc(clock);
                printf("[predict_stuff] Loaded predictions in %lfs\n", elapsed_time(clock));
                printf("[predict_stuff] Computing differences\n");
                FILE *diff_csv = fopen(diff_out, "w");
                fprintf(diff_csv, "diff\n");
                Vector_print_head_f(f_hat, 10);

                // Matrix_d *residuals = Matrix_new_d(1, Matrix_size_f(f_hat));

                // double sum_res = 0;
                double res_sq = 0;
                double abs_res = 0;


                for (int i = 0; i < args->data->N; i++) {
                    double res = f_hat->data[i] - (float) args->data->r[i];
                    // sum_res += res;
                    res_sq += res * res;
                    abs_res += fabs(res);

                    // printf("f_hat->data[i]: %f, data->r[i]: %u, residual: %f\n", f_hat->data[i], args->data->r[i], res);
                    // *vecacc_f(residuals, i) = res;
                    fprintf(diff_csv, "%f\n", res);
                }

                double mean_sq  = res_sq / args->data->N;
                double mean_abs = abs_res / args->data->N;

                fclose(diff_csv);
                
                printf("MSE: %lf\n", mean_sq);
                printf("AAE: %lf\n", mean_abs);

                Matrix_free_f(f_hat);

            // } else {

                // printf("[main] Differences already computed, stored in %s\n", diff_out);
                // fclose(test_open);
            // } // if (test_open == NULL)

            // With the differences computed, print a few statistics about the residuals

            

        } // if (diff)

        free(clock);
    } // if (predict)
}

void predict_stuff_mpi(const args_t *args, const SphericalModel *model, const Precomp *precomp) {

    // compute the predicted values and store them in a binary file
    // fhat_<size_dataset>_<lmodel>.bin
    int this_rank = args->rank;
    int world_size = args->world_size;


    char predictions_binary[100] = {0};
    sprintf(predictions_binary, "fhat_%s_%lu.bin", args->size_dataset, model->lmax);


    // Matrix_f *f_hat = NULL;
    Clock *clock = Clock_new();

    FILE *test_open = fopen(predictions_binary, "rb");
    if (test_open == NULL || args->recompute) {

        /**========================================================================
         *!                           Compute Prediction
            *========================================================================**/
        MPI_ONCE (
            printf("[predict_stuff_mpi] Computing altitude predictions...\n");
        )

        MPI_Barrier(MPI_COMM_WORLD);
        Clock_tic(clock);
        Matrix_f *f_hat = compute_prediction_mpi(model, precomp, args->data, args->size_dataset, world_size, this_rank);

        MPI_ONCE(

            printf("[predict_stuff_mpi] size f_hat: %lu\n", Matrix_size_f(f_hat));
        )

        Clock_toc(clock);
        if (this_rank == 0) {
            printf("[main] Time to compute prediction: %lfs\n", elapsed_time(clock));
            save_prediction(f_hat, predictions_binary);
            // Vector_print_head_f(f_hat, 10);
        }

    } else {

        MPI_ONCE (
            printf("[main] Predictions already computed\n");
        )
        fclose(test_open);
    }

    MPI_ONCE(
        if (args->print_args)
            head_data(args->data);
    )

    // We can't reach this point before having the predictions already computed
    if (args->diff) {

        char diff_out[100] = {0};
        sprintf(diff_out, "diff_%s_%d.csv", args->size_dataset, args->lmodel);

        // test_open = fopen(diff_out, "r");
        // if (test_open == NULL) {

            Clock_tic(clock);
            Matrix_f *f_hat = load_prediction(predictions_binary, args->data->N);
            Clock_toc(clock);

            FILE *diff_csv = NULL;
            MPI_ONCE(

                printf("[main] Loaded predictions in %lfs\n", elapsed_time(clock));
                printf("[main] Computing differences\n");
                diff_csv = fopen(diff_out, "w");
                double res_sq = 0;
                double abs_res = 0;


                for (int i = 0; i < args->data->N; i++) {
                    double res = f_hat->data[i] - (float) args->data->r[i];
                    // sum_res += res;
                    res_sq += res * res;
                    abs_res += fabs(res);

                    // printf("f_hat->data[i]: %f, data->r[i]: %u, residual: %f\n", f_hat->data[i], args->data->r[i], res);
                    // *vecacc_f(residuals, i) = res;
                    fprintf(diff_csv, "%f\n", res);
                }

                double mean_sq  = res_sq / args->data->N;
                double mean_abs = abs_res / args->data->N;

                fclose(diff_csv);
                
                printf("MSE: %lf\n", mean_sq);
                printf("AAE: %lf\n", mean_abs);

                Matrix_free_f(f_hat);
            )


        // } else {

            // printf("[main] Differences already computed, stored in %s\n", diff_out);
            // fclose(test_open);
        // } // if (test_open == NULL)
    } // if (diff)
}
