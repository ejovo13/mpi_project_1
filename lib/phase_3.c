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
