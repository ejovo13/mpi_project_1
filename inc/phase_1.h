#ifndef GEODESY_PHASE_1_H
#define GEODESY_PHASE_1_H

/**========================================================================
 * ?                          phase_1.h
 * @brief   : Compute Binary files to store the computations of Spherical
 *            harmonics for all i.
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-17
 *========================================================================**/
#include <stdio.h>

#include "constants.h"
#include "data.h"
#include "harmonics.h"
#include "learn.h"
#include "model.h"

// Compute the values [P00 P10 P11 .. Plm](cos \th) for all discrete values of theta
// and then store the results in a binary file whose first (lmax + 1)(lmax + 2)/2 * 8 bytes correspond to the 
// evaluation of all Plm at the first \th value.
/**
 * @brief Create a binary plm object
 *
 * @param lmax    The highest degree to precompute the values for.
 * @param th 
 */
void write_binary_plm(int lmax, const Matrix_d *theta, const char *binary_out);

// Return the entire P_lm_th matrix from a binary file 
Matrix_d *read_binary_plm(int lmax, int n_theta, const char *binary_filename);

// Return a matrix containing only the P_lm corresponding to l = L
// Return value is a n_theta x (L + 1) matrix
Matrix_d *read_binary_plm_l(int lmax, int l, int n_theta, const char *binary_filename);

// Need a RANGE function that returns a precomputed set of L's between a certain range

/**========================================================================
 *!               More Robust Precomputed Values
 *========================================================================**/
// A precomp structure stores and retrieves values that have been precomputed by 
// write_binary_plm and also stores precomputation of values like sin(m\phi), cos(m\phi),
// and sin(\theta)
typedef struct precomp_t {

    size_t L0; // Initial L value stored in Plm_th
    size_t LF; // Terminal L value stored in Plm_th
    Matrix_d *Plm_th;
    Matrix_d *sinth;  // sin(\th)     1 x n_theta
    Matrix_d *sinmph; // sin(m\phi)   m x n_phi
    Matrix_d *cosmph; // cos(m\phi)   m x n_phi
 
} Precomp;


Precomp *newPrecomp(int L0, int LF, int Lmax, const data_iso *data, const char *plm_bin);

void freePrecomp(Precomp *precomp);

/**========================================================================
 *!                           SphericalModel + Precomp + Data
 *========================================================================**/

// Compute the Coefficients Clm and Slm up until l = lmax, and store them in a newly 
// allocated SphericalModel
SphericalModel *newSphericalModel(int lmax, const data_iso *data, const Precomp *precomp);

double modelComputeCSlmPrecomp(SphericalModel *model, const data_iso *data, const Precomp *precomp);

// Write the coefficients of a spherical model to CSV
// Use the prefix to indicate what data set the spherical model was trained on
// pass "" as type if you don't want to specify
void SphericalModelToTXT(const SphericalModel *model, const char *type);

// Write the coefficients of a spherical model to a binary file
// void SphericalModelToBIN



#endif //GEODESY_PHASE_1_H