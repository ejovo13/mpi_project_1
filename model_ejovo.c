
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "harmonics.h"
#include "gmp.h"

#define LMAX 0

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

int main() {


    // Start toying with the small data set
    // data_points_iso *data = load_data_points_iso("ETOPO1_small.csv", 64800);

    const double th_0 = 0;
    const double th_f = PI;
    const double ph_0 = 0;
    const double ph_f = TWO_PI;
    const int t = 180;
    const int p = 360;

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

    data_iso *data = load_iso("ETOPO1_small.csv", 180, 360);

    // print ph
    // printf("{ ");
    // for (int i = 0; i < p; i++) {
    //     printf("%7.3lf ", data->ph[i]);
    // }
    // printf("}\n");

    const int nrows = (LMAX + 1) * (LMAX + 2) / 2;

    //
    double *P = (double *) malloc(sizeof(*P) * nrows * t);

	// for (int i = 0; i < data->npoint; i++) {
    // printf("{ ");
	for (int i = 0; i < t; i++) { // store this data in a more compact matrix

        // double th = th_0 + d_th * (i - 1);
        // printf("%7.3lf ", th[i]);
        // if (i % (i / 10) == 0) 
            // printf("Computing P values for i: %d\n", i);
        // P is a big matrix who stores the information as (p00(th0) p10(th0) )
		computeP(&model, &P[i * nrows], cos(data->th[i])); // should be mathematically equal to P(sin(phi_p)), angle named
                                               // phi in the project
        // P[PT(l, m)]

        // multiply by sin(phi) 
        for (int j = 0; j < nrows; j++) {
            P[i * nrows + j] *= sin(data->th[i]);
        } 


    }

    double c_integral = 0;
    double s_integral = 0;

    for (int l = 0; l <= LMAX; l++) {

        for (int m = 0; m <= l; m++) {

            // Now that I have the Associated legendre functions and the data efficiently loaded, let's write
            // the code to approximate the integrals
            // for (int j = 0; j < data->p; j++) {
            for (int j = 0; j < data->p; j++) {

                double ph_j = data->ph[j];

                for (int i = 0; i < data->t; i++) {

                    c_integral += data->r[j * data->p + i] * P[data->t * i + PT(l, m)] * cos(m * ph_j);
                    s_integral += data->r[j * data->p + i] * P[data->t * i + PT(l, m)] * sin(m * ph_j);
                }
            }

            c_integral *= (d_ph * d_th);
            s_integral *= (d_ph * d_th);

            // Now need to normalize with the factorial expression
        }
    }

    printf("Computing factorial:\n");

    mpz_t x, y;
    int f = 100;

    mpz_init_set_ui(x, f);
    mpz_init(y);
    mpz_factorial(y, x);

    gmp_printf("%d! = %.Zd\n", f, y);

    // printf("}\n");

    // print the matrix (at least the first 10 elements)
    // for (int i = 0; i < nrows; i++) {
    //     printf("| ");
    //     for (int j = 0; j < t; j++) {
    //         printf("%7.3lf ", P[j * nrows + i]);
    //     }
    //     printf("|\n");
    // }

    write_iso(data, "iso2.csv");

}