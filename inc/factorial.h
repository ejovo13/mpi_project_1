#ifndef GEODESY_FACTORIAL_H
#define GEODESY_FACTORIAL_H


#include <stdlib.h>
#include <stdio.h>

#include "gmp.h"


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

#endif // GEODESY_FACTORIAL_H