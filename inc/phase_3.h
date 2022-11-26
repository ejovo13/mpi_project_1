#ifndef GEODESY_PHASE_3
#define GEODESY_PHASE_3

/**========================================================================
 * ?                          phase_3.h
 * @brief   : Optimize validation routines by utilizing precomputed
 *            values and loading data sets from binary.
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-19
 *========================================================================**/

#include "ejovo.h"

#include "phase_1.h"
#include "phase_2.h"
#include "cli.h"

// A prediction uses floats (to reduce size by a factor of 2) and is simple a length N row 
// vector that is computed from a model, precomputed values, and a data set. A prediction 
// will be able to be stored efficiently in binary so that we can effectively compute the 
// MSE associated 
Matrix_f *compute_prediction(const SphericalModel *model, const Precomp *precomp, const data_iso* data);

Matrix_f *compute_prediction_omp(const SphericalModel *model, const Precomp *precomp, const data_iso* data);

Matrix_f *compute_prediction_mpi(const SphericalModel *model, const Precomp *precomp, const data_iso* data, const char *dataset_size, int world_size, int this_rank);

float compute_prediction_point(const SphericalModel *model, const Precomp *precomp, const data_iso *data, int i);

// Save a prediction to a binary file for fast computation of MSE
void save_prediction(const Matrix_f *f_hat, const char *binary_out);

// Load a prediction of N points stored in the binary file binary_in
Matrix_f *load_prediction(const char *binary_in, int N);

void predict_stuff(const args_t *args, const SphericalModel *model, const Precomp *precomp);

void predict_stuff_mpi(const args_t *args, const SphericalModel *model, const Precomp *precomp);

#endif // GEODESY_PHASE_3