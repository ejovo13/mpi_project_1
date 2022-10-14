
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ejovo.h"
#include "harmonics.h"
#include "gmp.h"

#define LMAX 50


// Given an lmax, compute the values of C_lm and S_lm and return the values
// in a Vector of the form
//
//  Coeff = [ C_lm S_lm ]
//
// where the length of Coeff is 2 * (lmax + 1) * (lmax + 2) / 2
iso_model *compute_model(int lmax, const char *filename, int n_theta, int n_phi) {

    printf("[compute_model] Constructing model { lmax: %d, t: %d, p: %d}\n", lmax, n_theta, n_phi);
    data_iso *data = load_iso(filename, n_theta, n_phi);
    iso_model *iso = newModel(data, lmax);
    writeModel(iso->model, data, "");    

return iso;

    // print out some of the recorded data
    // printf("phi: \n");
    // Matrix_print_d(data->ph);

    // printf("theta: \n");
    // Matrix_print_d(data->th);



    // // Calculate the MSE using sph_mod.
    // double mse = compute_mse(data, model);



    // printf("mse: %lf\n", mse);

    // // Let's try and adjust the parameters
    // double start_average_error = compute_average_error(data, model);
    // double start_mse = compute_mse(data, model);


    // for (int i = 0; i < 100; i++) {

    //     Matrix_d *grad = compute_gradient(data, model);

    //     // normalize the gradient

    //     adjust_parameters(data, grad, model, 0.000005);
    //     printf("mse: %lf\n", compute_mse(data, model));

    //     free(grad);
    // }

    // double end_average_error = compute_average_error(data, model); 
    // double end_mse = compute_mse(data, model); 

    // printf("err_0 := %lf -> err_f := %lf\n", start_average_error, end_average_error);
    // printf("mse_0 := %lf -> mse_f := %lf\n", start_mse, end_mse);

    // writeModel(model, data, "optimised");

}

// const int nrows = (LMAX + 1) * (LMAX + 2) / 2;

Matrix_d *modelPredictN(const iso_model *iso, int __n);

// conversion functions from iso to project
// project to wolfram
double phip_to_thew(double phi) {
    return (PI / 2) - phi;
}

double lamp_to_phiw(double lambda) {
    if (lambda < 0) return lambda + TWO_PI;
    return lambda;
}

double thew_to_phip(double theta) {
    return (PI / 2) - theta;
}

double phiw_to_lamp(double phiw) {
    if (phiw >= PI) return phiw - TWO_PI;
    // return TWO_PI - phiw;
    return phiw;
}


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

    // iso_model *iso = compute_model(5, filename_s, 4, 5);
    iso_model *iso = compute_model(4, filename_u, t_u, p_u);
    // write_iso(iso->data, "iso.csv");

    // int n = 10;

    // print phi
    Matrix_print_d(iso->data->ph);
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
    Matrix_print_d(iso->data->r);

    // Matrix_print_d(iso->model->C_lm);
    // Matrix_print_d(iso->model->S_lm);

    printf("MSE(iso) : %lf\n", compute_mse(iso));
    printf("err(iso) : %lf\n", compute_average_error(iso));

    printf("\t========================f_hat ====================\t\n"); Matrix_print_d(f_hat);
    
    // Coefficients of Clm and Slm
    // printf("\t======================= C_lm := ===================\t\n");
    // Matrix_print_d(iso->model->C_lm);
    // printf("\t======================= S_lm := ===================\t\n");
    // Matrix_print_d(iso->model->S_lm);

    printf("P_lm_th: ");
    Matrix_print_d(iso->model->P_lm_th);
    // compute_model(5, filename_s, 5, 10);
}

/**========================================================================
 *!                           Multiple precision functions
 *========================================================================**/
// unary function
void mpz_factorial(mpz_t rop, const mpz_t op) {
    // first create a grid of factorials that 
    // are the same precision as rop.
    // int prec = rop->_mp_prec;

    // create a grid of length rop
    // make sure that the value is under 1000
    // if (mpz_cmp_d(op, 1000) > 0) err(1, "Factorial is too large");
    if (mpz_cmp_ui(op, 1000) > 0) exit(2);

    // gmp_printf("%.Ff", op);
    
    mpz_t *facs = (mpz_t *) malloc(sizeof(*facs) * mpz_get_ui(op));

    for (int i = 0; mpz_cmp_si(op, i) > 0; i++) {
        // printf("i: %d\n", i);
        // mpf_t tmp;
        // mpf_init_set_ui(tmp, i);
        // gmp_printf(" * %Ff", tmp);
        mpz_init_set_ui(facs[i], i + 1);
        gmp_printf(" * %Zd", facs[i]);
    }
    printf("\n");

    mpz_t tmp;
    mpz_init_set_ui(tmp, 1);

    for (int i = 0; mpz_cmp_si(op, i) > 0; i++) {
        mpz_mul(tmp, tmp, facs[i]);
        // mpf_clear(facs[i]);
    }

    mpz_set(rop, tmp);

}

// compute (l + m)(l + m - 1) ... (l)(l - 1)(l - 2) ... (l - m + 1)
// rop is assumed to have been initialized (mpz_init)
void mpz_compute_denom(mpz_t rop, int l, int m) {

    // init the result to 1
    mpz_set_ui(rop, 1);
    mpz_t temp;

    mpz_init(temp);
        
    // perform 2m calculations
    for (int i = 0; i < 2 * m; i++) {
        mpz_set_ui(temp, (l + m - i));
        mpz_mul(rop, rop, temp);
    }

    mpz_clear(temp);
}

// compute the rational value (2l + 1)(l - m)! / (l + m)!
// and return a double value multiplied by (-1 / 2 * pi)
double factorial_coeff(int l, int m) {

    // create mpz (l - m)
    mpz_t den, num;
    mpz_init(den);

    mpz_init_set_ui(num, 2 * l + 1);
    mpz_compute_denom(den, l, m);

    // create a rational that is (2l + 1) / den
    mpq_t rat;

    mpq_set_num(rat, num);
    mpq_set_den(rat, den);

    double out = mpq_get_d(rat);

    if (out < 1e-100) return 0;

    // mpz_clear(den);
    // mpz_clear(num);
    // mpq_clear(rat);

    return out * (-1.0 / TWO_PI);
}

// Return a newly allocated matrix containing [c_l^m and s_li^m]
// O(lmax^2 * data->N) in time
// the return is a (2 * data->N) x (lmax^2 / 2)
//
// [[c00 c10 c11 ... clm](1)
//  [c00 c10 c11 ... clm](2)
//  [ .      .        . ](.)
//  [ .       .       . ](.)
//  [c00       .     clm](i)
//  [s00 s10 s11 ... clm](1)
//  [s00 s10 s11 ... clm](2)
//  [ .      .        . ](.)
//  [ .       .       . ](.)
//  [ .        .      . ](.)
//  [s00             slm](i)]
//
Matrix_d *compute_mse_coeff(const data_iso *data, const spherical_model *model) {

    // const int lmax = Matrix_size_d(P_lm_th);
    // const int 
    // const int ll = (lmax + 1) * (lmax + 2) / 2;
    const int ll = model->ll;
    const int lmax = model->lmax;



    // double *CS = (double *) malloc(sizeof(*CS) * 2 * ll * data->N); // start off using full sized array
    Matrix_d *cs_mat = Matrix_new_d(data->N * 2, ll);
    double *C = cs_mat->data;
    // double *S = matacc_d(cs_mat, 1, 0); // halfway point

    int index = 0;;
    double cosmph = 0;
    double sinmph = 0;
    const int N = data->N;

    // for (int i = 0; i < data->N; i++) {
    for (int i = 0; i < data->N; i++) {

        for (int l = 0; l <= lmax; l++) {
            // th changes every single cycle
            for (int m = 0; m <= l; m++) {

                int j = PT(l, m);
                 
                int i_ph = i / data->t; // ph has floor division cycles
                int i_th = i % data->t;
                
                cosmph = cos(m * vecat_d(data->ph, i_ph));
                sinmph = sin(m * vecat_d(data->ph, i_ph)); 

                *matacc_d(cs_mat, i, j)     = matat_d(model->P_lm_th, i_ph, i_th) * cosmph;
                *matacc_d(cs_mat, i + N, j) = matat_d(model->P_lm_th, i_ph, i_th) * sinmph;

                // printf("(%d) Computed mse coefficients c(%d, %d)\n", i, l, m);
            }
        }
        // printf("(i = %d) clm_i = ", i);
        // MatIter_print_d(Matrix_row_begin_d(cs_mat, i), Matrix_row_end_d(cs_mat, i));
        
    }

    return cs_mat;
}

Matrix_d *compute_mse_coeff_clm(const data_iso *data, const spherical_model* model) {

    Matrix_d *clm = Matrix_new_d(data->N, model->ll);
    double cosmph = 0;

    printf("lmax: %lu\n", model->lmax);

    for (int i = 0; i < data->N; i++) {
        for (int l = 0; l <= model->lmax; l++) {
            // th changes every single cycle
            for (int m = 0; m <= l; m++) {

                int j = PT(l, m);
                 
                int i_ph = i / data->t; // ph has floor division cycles
                int i_th = i % data->t;

                // printf("j: %d\n", j);
                
                // cosmph = cos(m * data->ph[i_ph]); 
                cosmph = cos(m * vecat_d(data->ph, i_ph)); 
                *matacc_d(clm, i, j)     = matat_d(model->P_lm_th, i_th, j) * cosmph;
            }
        }
    }

    return clm;
}

Matrix_d *compute_mse_coeff_slm(const data_iso *data, const spherical_model* model) {

    Matrix_d *slm = Matrix_new_d(data->N, model->ll);
    double sinmph = 0;

    for (int i = 0; i < data->N; i++) {
        for (int l = 0; l <= model->lmax; l++) {
            // th changes every single cycle
            for (int m = 0; m <= l; m++) {

                int j = PT(l, m);
                 
                int i_ph = i / data->t; // ph has floor division cycles
                int i_th = i % data->t;
                
                // sinmph = cos(m * data->ph[i_ph]); 
                sinmph = sin(m * vecat_d(data->ph, i_ph)); 
                *matacc_d(slm, i, j)     = matat_d(model->P_lm_th, i_th, j) * sinmph;
            }
        }
    }

    return slm;
}


// O(lmax^2)
double compute_pcs_i(const data_iso *data, const spherical_model *model, const Matrix_d *clm, const Matrix_d *slm, int i) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    // I already have the values of c_li^m and s_li^m, just need to
    // combine them with the parameters of my model
    double sum = 0;

    for (int l = 0; l <= model->lmax; l++) {
        for (int m = 0; m <= l; m++) {

            int index = PT(l, m);

            sum += model_C(model, l, m) * matat_d(clm, i, index);
            sum += model_S(model, l, m) * matat_d(slm, i, index);
        }
    }

    return sum;
}

Matrix_d *compute_pcs(const data_iso *data, const spherical_model *model, const Matrix_d *clm, const Matrix_d *slm) {

    // Pcs will simply have length data->N
    // double * pcs = (double *) malloc(sizeof(*pcs) * data->N);
    Matrix_d *pcs_mat = Matrix_new_d(1, data->N);

    // for i in {1..N} compute PCS_i 
    for (int i = 0; i < data->N; i++) {
        pcs_mat->data[i] = compute_pcs_i(data, model, clm, slm, i);
    }

    return pcs_mat;
}

double compute_gradient_clm(const data_iso *data, const spherical_model *model, const Matrix_d *clm, const Matrix_d *pcs, int l, int m) {

    double sum = 0;

    for (int i = 0; i < data->N; i++) {

        double clmi = matat_d(clm, i, PT(l, m));
        sum += 2 * clmi * pcs->data[i] - 2 * clmi * vecat_d(data->r, i);
    }

    return sum;
}

double compute_gradient_slm(const data_iso *data, const spherical_model *model, const Matrix_d *slm, const Matrix_d *pcs, int l, int m) {


    double sum = 0;

    for (int i = 0; i < data->N; i++) {

        double slmi = matat_d(slm, i, PT(l, m));
        // sum += 2 * slmi * pcs->data[i] - 2 * slmi * data->r[i];
        sum += 2 * slmi * pcs->data[i] - 2 * slmi * vecat_d(data->r, i);
    }

    return sum;
}

// Compute gradient in the form 
// [dMSE/dc00 dMSE/dc10 ... dMSE/dclm
//  dMSE/ds00 dMSE/ds10 ... dMSE/dslm]
// Assume model has C_lm and S_lm and P_lm_th
Matrix_d *compute_gradient(const data_iso *data, const spherical_model *model) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    const Matrix_d *clm  = compute_mse_coeff_clm(data, model); // 2 x ll matrix
    const Matrix_d *slm  = compute_mse_coeff_slm(data, model); // 2 x ll matrix
    const Matrix_d *pcs = compute_pcs(data, model, clm, slm); // sum(clm * Plmcos * cos)_i for i in 1..N

    // allocate space for output array

    // double *grad = (double *) malloc(sizeof(*grad) * 2 * ll);
    Matrix_d *grad = Matrix_new_d(2, ll);

    for (int l = 0; l <= model->lmax; l++) {
        for (int m = 0; m <= l; m++) {

            grad->data[PT(l, m)] = compute_gradient_clm(data, model, clm, pcs, l, m);
            grad->data[PT(l, m) + ll] = compute_gradient_slm(data, model, slm, pcs, l, m);

        }
    }

    return grad;
}



// Assume that a model object has been fully initialized
double compute_mse(const iso_model *iso) {

    const spherical_model *model = iso->model;
    const data_iso        *data  = iso->data;


    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    double sum = 0;

    // Compute MSE by getting a list of predictions
    const Matrix_d *predictions = modelPredictN(iso, iso->data->N);

    // now compute the difference between the two
    for (int i = 0; i < data->N; i++) {
        double diff = vecat_d(predictions, i) - vecat_d(data->r, i);
        sum += diff * diff;
    }

    return sum / data->N;
}

double compute_average_error(const iso_model *iso) {

    const spherical_model *model = iso->model;
    const data_iso        *data  = iso->data;


    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    double sum = 0;

    // Compute MSE by getting a list of predictions
    const Matrix_d *predictions = modelPredictN(iso, iso->data->N);

    // now compute the difference between the two
    for (int i = 0; i < data->N; i++) {
        sum += vecat_d(predictions, i) - vecat_d(data->r, i);
    }

    return sum / data->N;
}

// Given a gradient matrix and a model, tweak the values of 
// C_lm to improve the score
// gradient has size (N * 2) x ((lmax + 1)(lmax + 2)/2)
void adjust_parameters(const data_iso *data, const Matrix_d *grad, spherical_model *model, double alpha) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    // create a scaled version of the gradient vector (not necessary, but since the gradient vector should 
    // be reasonably small, we can perform this operation with very little cost)
    Matrix_d *grad_scaled = Matrix_clone_d(grad); // FAST memory copy
    Matrix_mult_scalar_d(grad_scaled, alpha); // Just multiply the first row
    // Matrix_mult_k

    // loop through the gradient matrix and adjust the model

    // matadd_d(model->)
    MatIter_row_add_row_d(Matrix_row_begin_d(model->C_lm, 0), 
                          Matrix_row_end_d(model->C_lm, 0), 
                          Matrix_row_begin_d(grad_scaled, 0));

    MatIter_row_add_row_d(Matrix_row_begin_d(model->S_lm, 0),
                          Matrix_row_end_d(model->S_lm, 0),
                          Matrix_row_begin_d(grad_scaled, 1));
   
    // for (int l = 0; l <= model->lmax; l++) {
    //     for (int m = 0; m <= l; m++) {

    //         // printf ("C(%d, %d): %lf\n", l, m, model->C_lm[PT(l, m)]);
    //        model_C_plus(model, l, m, -alpha * matat_d(grad, 0, PT(l, m))) -alpha * matat_d(grad, 0, PT(l, m));
    //        model_S(model, l, m) += -alpha * matat_d(grad, 1, PT(l, m));
    //     }
    // }
}