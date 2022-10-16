#ifndef GEODESY_LEARN_H
#define GEODESY_LEARN_H

/**========================================================================
 * ?                          learn.h
 * @brief   : Functions dealing with learning, prediction, and analysis
 *            of spherical models 
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-15
 *========================================================================**/
#include "constants.h"
#include "harmonics.h"
#include "data.h"
#include "model.h"
#include "ejovo.h"

/**========================================================================
 *!                           Prediction
 *========================================================================**/
double modelPredict(const iso_model *iso, double theta, double phi);

double modelPredictDataPoint(const iso_model *iso, int i);

Matrix_d *modelPredictDataPoints(const iso_model *iso, Matrix_i *indices);

Matrix_d *modelPredictData(const iso_model *iso);

Matrix_d *modelPredictN(const iso_model *iso, int __n);

void update_prediction(iso_model *iso);

/**========================================================================
 *!                           Learning
 *========================================================================**/
iso_trained *train_model(iso_model *iso, int nrunsdouble, double alpha_zero);

iso_trained *train_model_both(iso_model *iso, int nruns, double alpha_zero);

iso_trained *train_model_stochastic(iso_model *iso, int nruns, int __ndraws, double alpha_zero); 

// Compute gradient in the form 
// [dMSE/dc00 dMSE/dc10 ... dMSE/dclm
//  dMSE/ds00 dMSE/ds10 ... dMSE/dslm]
// Assume model has C_lm and S_lm and P_lm_th
Matrix_d *compute_gradient(const data_iso *data, const spherical_model *model);

Matrix_d *compute_stochastic_gradient(const iso_model *iso, int n);

/**========================================================================
 *!                           Analysis
 *========================================================================**/
// Assume that a model object has been fully initialized
double compute_mse(const iso_model *iso);

// What do I need to create a new model? A set of data and an lmax
double estimate_mse(const iso_model *iso, int n);

double compute_average_error(const iso_model *iso);

void print_predictions(const iso_model *iso);

void report_training(const iso_trained *trained, const char *prefix);

void print_coeff(const iso_model *iso);

// Given a gradient matrix and a model, tweak the values of 
// C_lm to improve the score
// gradient has size (N * 2) x ((lmax + 1)(lmax + 2)/2)
void adjust_parameters(const Matrix_d *grad, spherical_model *model, double alpha);

#endif // GEODESY_LEARN_H