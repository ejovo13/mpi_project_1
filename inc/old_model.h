#ifndef GEODESY_OLD_MODEL_H
#define GEODESY_OLD_MODEL_H

/**========================================================================
 * ?                          old_model.h
 * @brief   : Cursed model using least squares that runs in O(Nlmax^4) yikes
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-18
 *========================================================================**/

#include "harmonics.h"
#include "data.h"

/**************************** LINEAR ALGEBRA *********************************/

/*
 * Return the euclidean norm of x[0:n] using tricks for a greater precision
 */
double norm(int n, double const *x);

/*
 * Apply a real elementary reflector H to a real m-by-n matrix C. H is
 * represented in the form
 *
 *       H = I - tau * v * v**T
 *
 * where tau is a real scalar and v is a real vector.
 *
 * C is a 2D array of dimension (ldc, n).  On exit, C is overwritten with H*C.
 * It is required that ldc >= m.
 */
void multiply_householder(int m, int n, double *v, double tau, double *c, int ldc);

/*
 * Compute a QR factorization of a real m-by-n matrix A (with m >= n).
 *
 * A = Q * ( R ),       where:        Q is a m-by-n orthogonal matrix
 *         ( 0 )                      R is an upper-triangular n-by-n matrix
 *                                    0 is a (m-n)-by-n zero matrix
 *
 * A is a 2D array of dimension (m, n)
 * On exit, the elements on and above the diagonal contain R; the elements below
 * the diagonal, with the array tau, represent the orthogonal matrix Q.
 *
 * Q is represented as a product of n elementary reflectors
 *
 *     Q = H(1) * H(2) * ... * H(n).
 *
 *  Each H(i) has the form
 *
 *     H(i) = I - tau[i] * v * v**T
 *
 * where tau[i] is a real scalar, and v is a real vector with v[0:i-1] = 0 and 
 * v[i] = 1; v[i+1:m] is stored on exit in A[i+1:m, i].
 */
void QR_factorize(int m, int n, double * A, double * tau);

/*
 * Overwrite vector c with transpose(Q) * c where Q is a
 * real m-by-m orthogonal matrix defined as the product of k elementary
 * reflectors
 *
 *       Q = H(1) * H(2) ... H(k)
 * 
 * A is a 2D array of dimension (m, k), which contains a QR factorisation
 * computed by QR_factorize().  A is not modified.
 *
 * tau is an array of dimension k. tau[i] must contain the scalar factor of the
 * elementary reflector H(i), as returned by QR_factorize().  tau is read-only.
 *
 * c is a vector of dimension m.  On exit, c is overwritten by transpose(Q)*c.
 */
void multiply_Qt(int m, int k, double * A, double * tau, double * c);

/*
 * Solve the triangular linear system U*x == b
 *
 * U is a 2D array of dimension (ldu, n) with non-zero diagonal entries. Only
 * the upper-triangle is read by this function. b and x are n element vectors.
 * On exit, b is overwritten with x.
 */
void triangular_solve(int n, const double *U, int ldu, double *b);

/*
 * Solve the least-Squares Problem min || A*x - b || for overdetermined real
 * linear systems involving an m-by-n matrix A using a QR factorization of A.
 * It is assumed that A has full rank (and m >= n).
 *
 * A is a 2D array of dimension (m, n).  On exit, A is overwritten by the 
 * details of its QR factorization (cf. QR_factorize).
 *
 * b is a vector of size m.  On exit, b[0:n] contain the least squares solution 
 * vector; the residual sum of squares for the solution is given by the sum of 
 * squares of b[n:m].
 */
void linear_least_squares(int m, int n, double *A, double *b);

/*****************************************************************************/

// Compute the old model for a given lmax and number of points and then return how long it took.
// The clock starts immediately before computation of the least squares solution. Thus, when it 
// comes time to benchmark the new model, we should start the clock before the computation of Clm
// and Slm.
double old_model(int lmax, int npoint, const char *data_filename, const char *model_filename);


double old_model_quiet(int lmax, int npoint, const char *data_filename, const char *model_filename);

#endif // GEODESY_OLD_MODEL_H