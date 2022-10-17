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
 * @param dataset A string designating the ETOPO1 data set. Acceptable values are "small", "med", "high", "ultra"
 */
void write_binary_plm(int lmax, const Matrix_d *theta, const char *dataset);

// Return the entire P_lm_th matrix from a binary file 
Matrix_d *read_binary_plm(int lmax, int n_theta, const char *binary_filename);

// Return a matrix containing only the P_lm corresponding to l = L
Matrix_d *read_binary_plm_l(int lmax, int l, int n_theta, const char *binary_filename);


#endif //GEODESY_PHASE_1_H