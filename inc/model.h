#ifndef GEODESY_MODEL_H
#define GEODESY_MODEL_H
/**========================================================================
 * ?                          model.h
 * @brief   : Functions and structures dealing with the computational
 *            aspects of creating geodesic models
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-16
 *========================================================================**/
#include <math.h>

#include "data.h"

// Barebones model that contains ONLY the C_lm and S_lm 
// coefficients (plus bookkeeping fields)
typedef struct SphericalModel {

    size_t lmax;
    size_t ll;
    Matrix_d *C_lm;
    Matrix_d *S_lm;

} SphericalModel;

typedef struct spherical_model {
    size_t lmax;
    size_t ll;           // (lmax + 1)(lmax + 2) / 2
    Matrix_d *C_lm;      // O(ll) in space; 
    Matrix_d *S_lm;      // O(ll) in space;
    Matrix_d *P_lm_th;   // O(ll * t) in space
    Matrix_d *pcs;       // O(N) in space
    Matrix_d *clm;       // O(N * ll) in space
    Matrix_d *slm;       // O(N * ll) in space
} spherical_model;

// A model that has been solved with iso data
typedef struct iso_model {

    data_iso *data;                    // O(N)
    spherical_model * model;           // 
    struct spherical_harmonics *coeff;
    Matrix_d *f_hat; // prediction, not mandatory to have though

} iso_model;

typedef struct iso_trained {

    int count;
    iso_model *iso;
    Matrix_d *Clm_history; // ll x count matrix containing the training data
    Matrix_d *Slm_history; // ll x count matrix 
    Matrix_d *MSE;
    Matrix_d *dMSE;

} iso_trained;

// Simple structure that encodes the tuple (l, m)
typedef struct lm_pair {
    int l;
    int m;
} lm_pair;

static inline int PLM(int l, int m, int th, int nrows) {
    int index = PT(l, m);
    return th * nrows + index; 
}

/**
 * @brief Convert a 0-based index i into an (l, m) pair
 * 
 * Using a simple completing the square technique, solve for the value of l and then
 * use that information to compute the m index
 * 
 * @param i 
 * @return lm_pair 
 */
static inline lm_pair i_to_lm(int i) {

    double two_LL = (i + 1) * 2;
    double lmax_guess = sqrt(two_LL + 0.25) - 1.5;
    double l = ceil(lmax_guess);

    double LL_prev = l * (l + 1) / 2;
    double m = i - LL_prev;

    // Creationg of lm_pair
    lm_pair pair = {.l = l, .m = m};

    return pair;
}

// Print basic info about the model then print the first n coeff
static inline void Model_head(const SphericalModel *model) {
   printf("Model of degree %lu\n", model->lmax);
   Vector_print_head_d(model->C_lm, 10); 
   Vector_print_head_d(model->S_lm, 10); 
}

spherical_model *copy_spherical_model(const spherical_model *rhs);

iso_model *copy_iso (const iso_model *rhs);

// copy the current state of this iso, deep copying all of its components
iso_trained *copy_trained_model (const iso_trained *rhs);


// Assume that the model has already been initialized
/**========================================================================
 *!                           Accessor functions
 *========================================================================**/
static inline void model_C_plus(spherical_model *model, int l, int m, double addition) { model->C_lm->data[PT(l, m)] += addition; }
static inline void model_set_C(spherical_model *model, int l, int m, double val) { model->C_lm->data[PT(l, m)] = val; }
static inline double model_C(const spherical_model *model, int l, int m) { return model->C_lm->data[PT(l, m)]; }
static inline double model_S(const spherical_model *model, int l, int m) { return model->S_lm->data[PT(l, m)]; }
static inline double model_P(spherical_model *model, int l, int m, int thi) { return matat_d(model->P_lm_th, thi, PT(l, m)); } 

/**========================================================================
 *!                  Spherical Coordinate Transformations
 *========================================================================**/
static inline double phip_to_thew(double phi) { return (PI / 2) - phi; }
static inline double lamp_to_phiw(double lambda) {  return lambda; }
static inline double thew_to_phip(double theta) { return (PI / 2) - theta; }
static inline double phiw_to_lamp(double phiw) { return phiw; }

/**========================================================================
 *!                           Indexing Functions
 *========================================================================**/
static inline int data_i_ph(const data_iso *data, int i) { return i / data->t; }
static inline int data_i_th(const data_iso *data, int i) { return i % data->t; }

iso_model *newModel(data_iso *data, int lmax);

void writeModel(const spherical_model *model, const data_iso *data, const char *prefix);
void writeModifiedModel(const spherical_model *model, const data_iso *data, const char *prefix);

iso_model *compute_model(int lmax, const char *filename, int npoint);
iso_model *compute_model_binary(int lmax, const char *binary_in, int npoints, bool __log);

void free_spherical_model(spherical_model *sph_model);
void modelFree(iso_model* iso);

iso_model *newModelTimed(data_iso *data, int lmax, double *time_taken);
double time_new_model(int lmax, int npoint, const char *data_filename);

//! O(lmax^2 * t)  
// Compute the values P_lm(cos \theta) for every value of theta in a data set.
// For example, in ETOPO1_small.csv we have a discretized grid (\phi, \theta) \in [-pi, pi] x [pi/2, -pi/2] 
// where n_theta = 180 and n_phi = 360.
// The coefficients of P_lm are stored in a t x ll matrix, where ll = (lmax + 1) * (lmax + 2) / 2 and t = number of theta
//
// [[P00 P10 P11 P20 P21 ... Plm](cos \theta_1)
//  [P00 P10 P11 P20 P21 ... Plm](cos \theta_2)
//  [ .        .              . ](     .      )
//  [ .         .             . ](     .      )
//  [ .          .            . ](     .      )
//  [P00 P10 P11 P20 P21 ... Plm](cos \theta_t)]
void modelComputePlm(iso_model *model, const data_iso *data);

//! O(lmax^2 * N)
// Compute the values of the Laplace series coefficients C_lm and S_lm via the midpoint rule and store
// them in a spherical_model.
// There are LL = (lmax + 1) * (lmax + 2) / 2 Clm coefficients and LL Slm coefficients,
// which are stored under the C_lm and S_lm fields of a spherical_model. Both of these coefficients are 
// stored as 1 x LL matrices of doubles (Matrix_d) 
// C_lm : [C00 C10 C11 C20 C21 C22 ... Clm]
// S_lm : [S00 S10 S11 S20 S21 S22 ... Slm]
// return the walltime it took to compute this value
double modelComputeCSlm(spherical_model *model, const data_iso *data);

/**========================================================================
 *!                          Coefficients functions
 *========================================================================**/
//! O(lmax^2 * N) in time
// Compute the coefficients clm(i) := Plm(cos th_i)cos(m ph_i) and
//                          slm(i) := Plm(cos th_i)sin(m ph_i)
// for i in {1 .. N}
// 
// The coefficients clm and slm are used to compute the gradient of the MSE function, 
// with respect to the set of model parameters (C00, C10, C11, ..., Clm, S00, S10, S11, ..., Slm). 
//
// clm(i) and slm(i) are computed for the ith data point, thus we compute 2 * LL * N 
// coefficients, stored in a (N * 2) x (LL) matrix
//
// [[c00 c10     ... clm](1)
//  [c00 c10     ... clm](2)
//  [ .      .        . ](.)
//  [ .       .       . ](.)
//  [ .        .      . ](.)
//  [c00 c10         clm](N)
//  [s00 s10     ... slm](1)
//  [s00 s10     ... slm](2)
//  [ .      .        . ](.)
//  [ .       .       . ](.)
//  [ .        .      . ](.)
//  [s00             slm](N)]
Matrix_d *compute_mse_coeff(const data_iso *data, const spherical_model *model);

//! O(lmax^2 * N) in time
// Compute the coefficients clm(i) := Plm(cos th_i)cos(m ph_i) 
// That are used to compute the gradient of the MSE function, with respect to the set
// of model parameters (C00, C10, C11, ..., Clm, S00, S10, S11, ..., Slm)
//
// Store the coefficients clm(i) in a N x LL matrix
// [[c00 c10 c11 ... clm](1)
//  [c00 c10 c11 ... clm](2)
//  [ .      .        . ](.)
//  [ .       .       . ](.)
//  [ .        .      . ](.)
//  [c00             clm](N)]
Matrix_d *compute_mse_coeff_clm(const data_iso *data, const spherical_model* model);

//! O(lmax^2 * N) in time
// Compute the coefficients slm(i) := Plm(cos th_i)sin(m ph_i) 
// That are used to compute the gradient of the MSE function, with respect to the set
// of model parameters (C00, C10, C11, ..., Clm, S00, S10, S11, ..., Slm)
//
// Store the coefficients clm(i) in a N x LL matrix
// [[s00 s10 s11 ... slm](1)
//  [s00 s10 s11 ... slm](2)
//  [ .      .        . ](.)
//  [ .       .       . ](.)
//  [ .        .      . ](.)
//  [s00             slm](N)
Matrix_d *compute_mse_coeff_slm(const data_iso *data, const spherical_model* model);

//! O(lmax^2) in time
// 
// Compute the ith value PCS(i) :=  (c00_i * C00_i + s00_i * C00_i) + 
//                                  (c10_i * C10_i + s10_i * S10_i) + 
//                                                .                 + 
//                                                .                 + 
//                                                .                 + 
//                                  (clm_i * Clm_i + slm_i * Slm_i) 
//
// Since the values clm(i) and slm(i) have been precomputed (in O(lmax^2 * N) time), we can access them in 
// constant time for this function and our time complexity is thus O(lmax^2)
double compute_pcs_i(const spherical_model *model, const Matrix_d *clm, const Matrix_d *slm, int i);

Matrix_d *compute_pcs_indices(const iso_model *iso, const Matrix_i *indices);

//! O(lmax^2 * N) in time
// 
// Compute the values PCS(i) :=  (c00_i * C00_i + s00_i * C00_i) + 
//                               (c10_i * C10_i + s10_i * S10_i) + 
//                                             .                 + 
//                                             .                 + 
//                                             .                 + 
//                               (clm_i * Clm_i + slm_i * Slm_i) 
// for i in {1 .. N}
//
// Since the values clm(i) and slm(i) have been precomputed (in O(lmax^2 * N) time), we can access them in 
// constant time for this function and our time complexity is thus O(lmax^2)
Matrix_d *compute_pcs(const data_iso *data, const spherical_model *model, const Matrix_d *clm, const Matrix_d *slm);

// Compute the gradient of the MSE with respect to (C00, C10, C11, ..., Clm)
double compute_gradient_clm(const data_iso *data, const Matrix_d *clm, const Matrix_d *pcs, int l, int m);

// Compute the gradient of the MSE with respect to (S00, S10, S11, ..., Slm)
double compute_gradient_slm(const data_iso *data, const Matrix_d *slm, const Matrix_d *pcs, int l, int m);

// Compute an estimate to the gradient of the MSE with respect (C00, C10, C11, ..., Clm) using the 
// points prescribed by the indices matrix
double compute_gradient_clm_points(const iso_model *iso, const Matrix_d *clm, const Matrix_d *pcs, int l, int m, const Matrix_i *indices);

// Compute an estimate to the gradient of the MSE with respect (S00, S10, S11, ..., Slm) using the 
// points prescribed by the indices matrix
double compute_gradient_slm_points(const iso_model *iso, const Matrix_d *slm, const Matrix_d *pcs, int l, int m, const Matrix_i *indices);

// Compute the mean squared error of a model
double compute_mse(const iso_model *model);






#endif // GEODESY_MODEL_H