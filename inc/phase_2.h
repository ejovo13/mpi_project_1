#ifndef GEODESY_PHASE_2_H
#define GEODESY_PHASE_2_H

/**========================================================================
 * ?                          phase_2.h
 * @brief   : Assuming that an appropriate binary file has been computed, 
 *            output the coefficients Clm and Slm for a SINGLE value of L. 
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-17
 *========================================================================**/
#include <stdint.h>

#include "data.h"
#include "phase_1.h"

int16_t *read_binary_dataset(int n_points, const char *binary_file_in);

void reduce_csv_to_binary(int n_points, const char *csv_in, const char *binary_out);

// Write binary file containing the altitude values to a csv file for visual inspection
void write_binary_z_to_csv(int n_points, const int16_t *z, const char *output_csv);

void binary_dataset_to_csv(int n_points, const char *binary_in, const char *output_csv);

// Compute a full set of coefficients CLm and SLm, then store them to two binary files 
void write_binary_cslm(int L, int lmax, int n_theta, const data_iso *data, const char *binary_file_in);



#endif // GEODESY_PHASE_2_H