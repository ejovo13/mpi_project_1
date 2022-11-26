#include "geodesy.h"
#include "assert.h"

int main() {

    // Here we are going to test my i_to_lm function

    int lmax = 500;

    int LL = (lmax + 1) * (lmax + 2) / 2;

    Matrix_i *L = Matrix_new_i(1, LL);
    Matrix_i *M = Matrix_new_i(1, LL);

    // Properly fill the matrices
    int i = 0;

    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            
            vecset_i(L, i, l);
            vecset_i(M, i, m);

            i++;
        }
    }

    // Matrix_i *m_pair = Matrix_new_i(1, LL);
    // Now test to make sure the values are properly returned by i_to_lm

    for (int i = 0; i < LL; i++) {

        lm_pair pair = i_to_lm(i);
        assert(pair.l == vecat_i(L, i)); // The l's have been properly computed
        assert(pair.m == vecat_i(M, i));
    }

    Matrix_free_i(L);
    Matrix_free_i(M);

    fprintf(stderr, "Assertions complete\n");

    return 0;
}