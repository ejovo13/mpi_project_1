
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "harmonics.h"
#include "gmp.h"

#define LMAX 50

// compute the coefficient using the gmp multi-precision library
// mpf_t compute_scale(int l, int m) {







// }

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

const int nrows = (LMAX + 1) * (LMAX + 2) / 2;

static inline int PLM(int l, int m, int th) {
    int index = PT(l, m);
    return th * nrows + index; 
}


int main(int argc, char *argv[]) {


    // Start toying with the small data set
    // data_points_iso *data = load_data_points_iso("ETOPO1_small.csv", 64800);
    // int LMAX = 0;
    // if (argc >= 2) LMAX = 
    



    const double th_0 = 0;
    const double th_f = PI;
    const double ph_0 = 0;
    const double ph_f = TWO_PI;
    const int t = 180;
    const int p = 360;

    // double C_lm = double 

    const double d_th = (th_f - th_0) / t;
    const double d_ph = (ph_f - ph_0) / p;

    printf("th: [%lf, %lf]\tph: [%lf, %lf]\n", th_0, th_f, ph_0, ph_f);

    // printf("My model\n");

    // Start off by loading in the ETOPO1_small.csv data set

    // data_points_iso *data = load_data_points_iso("ETOPO1_small.csv", 100);

    // print_npoints(data, 10);
    // write_npoints(data, 64800, "iso.csv");
    // printf("Data has %d points\n", data->npoint);

    struct spherical_harmonics model;
	setup_spherical_harmonics(LMAX, &model);

    data_iso *data = load_iso("ETOPO1_small.csv", t, p);

    // print ph
    // printf("{ ");
    // for (int i = 0; i < p; i++) {
    //     printf("%7.3lf ", data->ph[i]);
    // }
    // printf("}\n");

    // const int nrows = (LMAX + 1) * (LMAX + 2) / 2;



    //
    double *P_lm_th = (double *) malloc(sizeof(*P_lm_th) * nrows * t);
    double *C_lm    = (double *) malloc(sizeof(*C_lm) * nrows);
    double *S_lm    = (double *) malloc(sizeof(*S_lm) * nrows);

	// for (int i = 0; i < data->npoint; i++) {
    // printf("{ ");
	for (int i = 0; i < t; i++) { // store this data in a more compact matrix

        // double th = th_0 + d_th * (i - 1);
        // printf("%7.3lf ", th[i]);
        // if (i % (i / 10) == 0) 
            // printf("Computing P values for i: %d\n", i);
        // P is a big matrix who stores the information as (p00(th0) p10(th0) )
		computeP(&model, &P_lm_th[i * nrows], cos(data->th[i])); // should be mathematically equal to P(sin(phi_p)), angle named
                                               // phi in the project
        // P[PT(l, m)]

        // multiply by sin(phi) 
        for (int j = 0; j < nrows; j++) {
            P_lm_th[i * nrows + j] *= sin(data->th[i]);
        } 


    }

    double c_integral = 0;
    double s_integral = 0;
    int count = 0;

    for (int l = 0; l <= LMAX; l++) {

        for (int m = 0; m <= l; m++) {

            // Now that I have the Associated legendre functions and the data efficiently loaded, let's write
            // the code to approximate the integrals
            // for (int j = 0; j < data->p; j++) {
            // int index = PT(l, m);

            c_integral = 0;
            s_integral = 0;

            for (int j = 0; j < data->p; j++) {

                double ph_j = data->ph[j];
                double cos_mph = cos(m * ph_j);
                double sin_mph = sin(m * ph_j);

                for (int i = 0; i < data->t; i++) {

                    // printf("(%d) data->t: %d\n", i, data->t); 

                    // printf("data->t * i + index = %d * %d + %d = %d\n", data->t, i, index, data->t * i + index);
                    // if (j == 0) 
                        // printf("(%d) P_%d^%d(th_%d) = %lf\n", index, l, m, i, P_lm_th[PLM(l, m, i)]);

                    // printf("(%d) data: %lf\n", count, data->r[j * data->t + i]);
                    c_integral += data->r[j * data->t + i] * P_lm_th[PLM(l, m, i)] * cos_mph;
                    s_integral += data->r[j * data->t + i] * P_lm_th[PLM(l, m, i)] * sin_mph;

                    count ++;
                }
            }

            // printf("===========================\n");
            c_integral *= (d_ph * d_th) / (TWO_PI);
            s_integral *= (d_ph * d_th) / (TWO_PI);


            // Now need to normalize with the factorial expression
            // But I actually think that this shit is already normalized
            C_lm[PT(l, m)] = c_integral;
            S_lm[PT(l, m)] = s_integral;

            // let's print out the computed values
            printf("C(%d, %d) = %lf \t S(%d, %d) = %lf\n", l, m, c_integral, l, m, s_integral);
            // printf("C(%d, %d) = %lf\n", l, m, c_integral);
            // printf("S(%d, %d) = %lf\n", l, m, s_integral);
            // printf("C(%d, %d) = %lf\n", l, m, c_integral);
        }
    }

    // output C_lm and S_lm to a file
    const char filename[] = "ejovo50.txt";

    FILE *out = fopen(filename, "w");

    for (int l = 0; l <= LMAX; l++) {
        for (int m = 0; m <= l; m++) {
            fprintf(out, "%d\t%d\t%.15lf\t%.15lf\n", l, m, C_lm[PT(l, m)], S_lm[PT(l, m)]);
        }
    }

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

    printf("=============================\n");

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

    write_iso(data, "iso2.csv");

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