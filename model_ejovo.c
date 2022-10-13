
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "harmonics.h"
#include "gmp.h"

#define LMAX 50

// compute the coefficient using the gmp multi-precision library
// mpf_t compute_scale(int l, int m) {

#include "ejovo.h"




double compute_mse(const data_iso *data, const spherical_model *model);
double *compute_gradient(const data_iso *data, const spherical_model *model);
void adjust_parameters(const data_iso *data, const double *grad, spherical_model *model, double alpha);
double compute_average_error(const data_iso *data, const spherical_model *model);
// double compute_gradient_slm(const data_iso *data, const spherical_model *model, const double *s, const double *pcs, int l, int m);
// double compute_gradient_clm(const data_iso *data, const spherical_model *model, const double *c, const double *pcs, int l, int m);
// double *compute_pcs(const data_iso *data, const spherical_model *model, const double *cs);



// }

static inline int PLM(int l, int m, int th, int nrows) {
    int index = PT(l, m);
    return th * nrows + index; 
}

// Given an lmax, compute the values of C_lm and S_lm and return the values
// in a Vector of the form
//
//  Coeff = [ C_lm S_lm ]
//
// where the length of Coeff is 2 * (lmax + 1) * (lmax + 2) / 2
void compute_model(int lmax, const char *filename, int n_theta, int n_phi) {

    const double th_0 = 0;
    const double th_f = PI;
    const double ph_0 = 0;
    const double ph_f = TWO_PI;
    const int nrows = (lmax + 1) * (lmax + 2) / 2;

    // allocate the coeffient matrices
    double *CS_lm   = (double *) malloc(sizeof(*CS_lm) * 2 * nrows);
    double *P_lm_th = (double *) malloc(sizeof(*P_lm_th) * nrows * n_theta);

    double *C_lm = CS_lm;
    double *S_lm = &CS_lm[nrows]; 

    const double d_th = (th_f - th_0) / n_theta;
    const double d_ph = (ph_f - ph_0) / n_phi;

    // Load data
    data_iso *data = load_iso(filename, n_theta, n_phi);

    // // print the theta vector
    // printf("{");
    // for (int i = 0; i < n_theta - 1; i++) {
    //     printf("%lf, ", data->th[i]);
    // }
    // printf("%lf}\n", data->th[n_theta - 1]);

    // printf("==================================\n");

    // // print phi
    // printf("{");
    // for (int i = 0; i < n_phi - 1; i++) {
    //     printf("%lf, ", data->ph[i]);
    // }
    // printf("%lf}\n", data->ph[n_phi - 1]);



    // Initialize P_lm computation coefficients
    struct spherical_harmonics model;
	setup_spherical_harmonics(lmax, &model);

    printf("th: [%lf, %lf]\tph: [%lf, %lf]\n", th_0, th_f, ph_0, ph_f);

    // make sure data got loaded
    printf("Data: {.N = %d, .t = %d, .p = %d}", data->N, data->t, data->p);

	for (int i = 0; i < n_theta; i++) { // store this data in a more compact matrix

        // P is a big matrix who stores the information as (p00(th0) p10(th0) )
		computeP(&model, &P_lm_th[i * nrows], cos(data->th[i])); // should be mathematically equal to P(sin(phi_p)), angle named
		// computeP(&model, &P_lm_th[i * nrows], sin(data->th[i])); // should be mathematically equal to P(sin(phi_p)), angle named
                                               // phi in the project

        // multiply by sin(phi) 
        for (int j = 0; j < nrows; j++) {
            P_lm_th[i * nrows + j] *= sin(data->th[i]);
            // P_lm_th[i * nrows + j] *= cos(data->th[i]);
        } 
    }

    double c_integral = 0;
    double s_integral = 0;
    int count = 0;

    for (int l = 0; l <= lmax; l++) {

        for (int m = 0; m <= l; m++) {

            // Now that I have the Associated legendre functions and the data efficiently loaded, let's write
            // the code to approximate the integrals
            c_integral = 0;
            s_integral = 0;

            for (int j = 0; j < data->p; j++) {

                double ph_j = data->ph[j];
                double cos_mph = cos(m * ph_j);
                double sin_mph = sin(m * ph_j);

                for (int i = 0; i < data->t; i++) {

                    c_integral += data->r[j * data->t + i] * P_lm_th[PLM(l, m, i, nrows)] * cos_mph;
                    s_integral += data->r[j * data->t + i] * P_lm_th[PLM(l, m, i, nrows)] * sin_mph;

                    count ++;
                }
            }

            // printf("===========================\n");
            c_integral *= (d_ph * d_th) / (TWO_PI);
            // c_integral *= (d_ph * d_th) / (PI);
            s_integral *= (d_ph * d_th) / (TWO_PI);
            // s_integral *= (d_ph * d_th) / (PI);


            // But I actually think that this shit is already normalized
            C_lm[PT(l, m)] = c_integral;
            S_lm[PT(l, m)] = s_integral;

        }
    }

    // output C_lm and S_lm to a file
    // const char filename[] = "ejovo50.txt";

    char file_out[100] = {0};

    //model_lmax_nth_nph.txt
    sprintf(file_out, "model_%d_%d_%d.txt", lmax, n_theta, n_phi);

    printf("\nWriting data to file: %s\n\n", file_out);

    FILE *out = fopen(file_out, "w");

    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            fprintf(out, "%d\t%d\t%.15lf\t%.15lf\n", l, m, C_lm[PT(l, m)], S_lm[PT(l, m)]);
        }
    }

    fclose(out);


    /**========================================================================
     *!                           Optimization 
     *========================================================================**/
    // Now that we've computed C_lm and S_lm, let's go ahead and store them in a model object.
    spherical_model sph_mod = {

        .lmax = model.lmax,
        .C_lm = C_lm,
        .S_lm = S_lm,
        .P_lm_th = P_lm_th

    };

    // Calculate the MSE using sph_mod.

    double mse = compute_mse(data, &sph_mod);

    printf("mse: %lf\n", mse);

    // Let's get wild and compute the gradient

    // double *grad = compute_gradient(data, &sph_mod);

    // what's the gradient look like?

    // for (int l = 0; l <= lmax; l++) {
    //     for (int m = 0; m <= l; m++) {
    //         printf("gradC(%d, %d) := %lf\n", l, m, grad[PT(l, m)]);
    //     }
    // }

    //

    // Let's try and adjust the parameters
    double start_average_error = compute_average_error(data, &sph_mod);
    double start_mse = compute_mse(data, &sph_mod);


    for (int i = 0; i < 100; i++) {

        double *grad = compute_gradient(data, &sph_mod);
        adjust_parameters(data, grad, &sph_mod, 0.000005);
        printf("mse: %lf\n", compute_mse(data, &sph_mod));

        free(grad);
    }

    double end_average_error = compute_average_error(data, &sph_mod); 
    double end_mse = compute_mse(data, &sph_mod); 

    printf("err_0 := %lf -> err_f := %lf\n", start_average_error, end_average_error);
    printf("mse_0 := %lf -> mse_f := %lf\n", start_mse, end_mse);

    sprintf(file_out, "adjusted_model_%d_%d_%d.txt", lmax, n_theta, n_phi);

    printf("\nWriting data to file: %s\n\n", file_out);

    out = fopen(file_out, "w");

    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            fprintf(out, "%d\t%d\t%.15lf\t%.15lf\n", l, m, sph_mod.C_lm[PT(l, m)], sph_mod.S_lm[PT(l, m)]);
        }
    }
    
    fclose(out);

}

// const int nrows = (LMAX + 1) * (LMAX + 2) / 2;




int main(int argc, char *argv[]) {


    // Start toying with the small data set
    // data_points_iso *data = load_data_points_iso("ETOPO1_small.csv", 64800);
    // int LMAX = 0;
    // if (argc >= 2) LMAX = 
    

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

    // Now create an interface for the function that computes 
    // int lmax = 10;

    // compute_model(lmax, filename_m, t_m, p_m); 

    // for (int l = 0; l < lmax; l++) {
        // compute_model(l, filename_s, t_s, p_s);
    // }
    // int l = 200;

    compute_model(5, filename_s, t_s, p_s);
    // compute_model(20, filename_m, t_m, p_m);

    // printf("Computing factorial:\n");

    // mpz_t x, y;
    // int f = 100;

    // mpz_init_set_ui(x, f);
    // mpz_init(y);
    // mpz_factorial(y, x);

    // gmp_printf("%d! = %.Zd\n", f, y);

    // printf("}\n");

    // // print the matrix (at least the first 10 elements)
    // for (int i = 0; i < nrows; i++) {
    //     printf("| ");
    //     for (int j = 0; j < t; j++) {
    //         printf("%7.3lf ", P_lm_th[j * nrows + i]);
    //         // printf("%7.3lf ", P_lm_th[]);
    //     }
    //     printf("|\n");
    // }

    // printf("=============================\n");

    // for (int l = 0; l <= LMAX; l++) {
    //     for (int m = 0; m <= l; m++) {

    //         printf("(%d, %d) ", l, m);
    //         printf("| ");
    //         for (int j = 0; j < t; j++) {
    //             printf("%7.3lf ", P_lm_th[PLM(l, m, j)]);
    //             // printf("%7.3lf ", P_lm_th[]);
    //         }
    //         printf("|\n");
    //     }
    // }

    // write_iso(data, "iso2.csv");

    // for a given l, let's compute the denominator
    // int l = 17;
    // mpz_t den;
    // mpz_init(den);

    // for (int m = 0; m <= l; m++) {

    //     mpz_compute_denom(den, l, m);
    //     gmp_printf("den(%d, %d) = %Zd\n", l, m, den);



    // }

    // Now just compute the factorial coefficient
    // for (int m = 0; m <= l; m++) {
    //     printf("Coeff(%d, %d): %lf\n", l, m, factorial_coeff(l, m));
    // }

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
double *compute_mse_coeff(const data_iso *data, const double *P_lm_th, int lmax) {

    const int ll = (lmax + 1) * (lmax + 2) / 2;

    double *CS = (double *) malloc(sizeof(*CS) * 2 * ll * data->N); // start off using full sized array
    double *C = CS;
    double *S = &CS[ll * data->N]; // halway point

    int index;
    double cosmph;
    double sinmph;

    // for (int i = 0; i < data->N; i++) {
    for (int i = 0; i < 300; i++) {
    // for (int i = 0; i < 20; i++) {


        for (int l = 0; l <= lmax; l++) {
            // th changes every single cycle
            for (int m = 0; m <= l; m++) {

                index = PT(l, m);
                 

                int i_ph = index / data->p; // ph has floor division cycles
                int i_th = index % data->t;
                
                cosmph = cos(m * data->ph[i_ph]); 
                sinmph = sin(m * data->ph[i_ph]); 

                C[i * ll + index] = P_lm_th[PLM(l, m, i_th, ll)] * cosmph;
                S[i * ll + index] = P_lm_th[PLM(l, m, i_th, ll)] * sinmph;

                // printf ("C(%d, %d, %d) := %lf\n", l, m, i, C[i * ll + index]); 
            }
        }
        // printf("(%d) Computed mse coefficients c(%d, %d)\n", i, l, m);
        // printf("(%d) Computed mse coefficients\n", i);
    }

    return CS;
}

// O(lmax^2)
double compute_pcs_i(const data_iso *data, const spherical_model *model, const double *cs, int i) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    // for i in {1..N} compute PCS_i 
    const double *c = cs;
    const double *s = &cs[ll * data->N]; // halfway point

 
    // I already have the values of c_li^m and s_li^m, just need to
    // combine them with the parameters of my model

    double sum = 0;

    for (int l = 0; l <= model->lmax; l++) {
        for (int m = 0; m <= l; m++) {

            int index = PT(l, m);

            sum += model->C_lm[index] * c[i * ll + index];
            sum += model->S_lm[index] * s[i * ll + index];

        }
    }

    return sum;
}

double *compute_pcs(const data_iso *data, const spherical_model *model, const double *cs) {

    // Pcs will simply have length data->N
    double * pcs = (double *) malloc(sizeof(*pcs) * data->N);

    for (int i = 0; i < data->N; i++) {
        pcs[i] = compute_pcs_i(data, model, cs, i);
    }

    return pcs;
}

double compute_gradient_clm(const data_iso *data, const spherical_model *model, const double *c, const double *pcs, int l, int m) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    double sum = 0;

    for (int i = 0; i < data->N; i++) {

        double clmi = c[PT(l, m) + i * ll];
        sum += 2 * clmi * pcs[i] - 2 * clmi * data->r[i];

    }

    return sum;
}

double compute_gradient_slm(const data_iso *data, const spherical_model *model, const double *s, const double *pcs, int l, int m) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    double sum = 0;

    for (int i = 0; i < data->N; i++) {

        double slmi = s[PT(l, m) + i * ll];
        sum += 2 * slmi * pcs[i] - 2 * slmi * data->r[i];

    }

    return sum;
}

// Compute gradient in the form 
// [dMSE/dc00 dMSE/dc10 ... dMSE/dclm
//  dMSE/ds00 dMSE/ds10 ... dMSE/dslm]
// Assume model has C_lm and S_lm and P_lm_th
double *compute_gradient(const data_iso *data, const spherical_model *model) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    const double *cs  = compute_mse_coeff(data, model->P_lm_th, model->lmax);
    const double *pcs = compute_pcs(data, model, cs);

    const double *c = cs;
    const double *s = &cs[ll * data->N];

    // allocate space for output array

    double *grad = (double *) malloc(sizeof(*grad) * 2 * ll);

    for (int l = 0; l <= model->lmax; l++) {
        for (int m = 0; m <= l; m++) {

            grad[PT(l, m)] = compute_gradient_clm(data, model, c, pcs, l, m);
            grad[PT(l, m) + ll] = compute_gradient_slm(data, model, s, pcs, l, m);

        }
    }

    return grad;

}

// Assume that a model object has been fully initialized
double compute_mse(const data_iso *data, const spherical_model *model) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    // compute pcs
    // double pcs;
    double *cs = compute_mse_coeff(data, model->P_lm_th, model->lmax);
    double *c = cs;
    double *s = &cs[ll * data->N];
    // double *pcs = compute_pcs(data, model, cs);

    // 

    // compute the sum from i to N



    double sum = 0;

    // for (int i = 0; i < data->N; i++) {
    for (int i = 0; i < data->N; i++) {

        double inner = 0;

        for (int l = 0; l <= model->lmax; l++) {
            for (int m = 0; m <= l; m++) {

                // printf("+ %lf\n", model->C_lm[PT(l, m)]);
                inner += model->C_lm[PT(l, m)] * c[PT(l, m)];
                inner += model->S_lm[PT(l, m)] * s[PT(l, m)];
            }
        }


        sum += (inner - data->r[i]) * (inner - data->r[i]);
        // sum += (inner - data->r[i]);

    }

    return sum / data->N;
}

double compute_average_error(const data_iso *data, const spherical_model *model) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    // compute pcs
    // double pcs;
    double *cs = compute_mse_coeff(data, model->P_lm_th, model->lmax);
    double *c = cs;
    double *s = &cs[ll * data->N];
    // double *pcs = compute_pcs(data, model, cs);

    // 

    // compute the sum from i to N



    double sum = 0;

    // for (int i = 0; i < data->N; i++) {
    for (int i = 0; i < data->N; i++) {

        double inner = 0;

        for (int l = 0; l <= model->lmax; l++) {
            for (int m = 0; m <= l; m++) {

                // printf("+ %lf\n", model->C_lm[PT(l, m)]);
                inner += model->C_lm[PT(l, m)] * c[PT(l, m)];
                inner += model->S_lm[PT(l, m)] * s[PT(l, m)];
            }
        }


        // sum += (inner - data->r[i]) * (inner - data->r[i]);
        sum += (inner - data->r[i]);

    }

    return sum / data->N;
}

// Given a gradient matrix and a model, tweak the values of 
// C_lm to improve the score
// gradient has size (N * 2) x ((lmax + 1)(lmax + 2)/2)
void adjust_parameters(const data_iso *data, const double *grad, spherical_model *model, double alpha) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    // loop through the gradient matrix and adjust the model
    const double *grad_c = grad;
    const double *grad_s = &grad[ll];

    for (int l = 0; l <= model->lmax; l++) {
        for (int m = 0; m <= l; m++) {

            // printf ("C(%d, %d): %lf\n", l, m, model->C_lm[PT(l, m)]);


           model->C_lm[PT(l, m)] += -alpha * grad_c[PT(l, m)];
           model->S_lm[PT(l, m)] += -alpha * grad_s[PT(l, m)];
        }
    }
}