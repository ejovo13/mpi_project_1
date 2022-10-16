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
#include "data.h"



typedef struct spherical_model {
    size_t lmax;
    size_t ll;   // (lmax + 1)(lmax + 2) / 2
    Matrix_d *C_lm; // coefficients of the MODEL
    Matrix_d *S_lm;
    Matrix_d *P_lm_th;
    Matrix_d *pcs;
    Matrix_d *clm;
    Matrix_d *slm; // small coefficients used during computation of MSE, gradient
} spherical_model;

// A model that has been solved with iso data
typedef struct iso_model {

    data_iso *data;
    spherical_model * model;
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


spherical_model *copy_spherical_model(const spherical_model *rhs);

iso_model *copy_iso (const iso_model *rhs);

// copy the current state of this iso, deep copying all of its components
iso_trained *copy_trained_model (const iso_trained *rhs);

static inline int PLM(int l, int m, int th, int nrows) {
    int index = PT(l, m);
    return th * nrows + index; 
}

// Assume that the model has already been initialized
static inline void model_C_plus(spherical_model *model, int l, int m, double addition) { model->C_lm->data[PT(l, m)] += addition; }
static inline void model_set_C(spherical_model *model, int l, int m, double val) { model->C_lm->data[PT(l, m)] = val; }
static inline double model_C(const spherical_model *model, int l, int m) { return model->C_lm->data[PT(l, m)]; }
static inline double model_S(const spherical_model *model, int l, int m) { return model->S_lm->data[PT(l, m)]; }
static inline double model_P(spherical_model *model, int l, int m, int thi) { return matat_d(model->P_lm_th, thi, PT(l, m)); } 

// conversion functions from iso to project
// project to wolfram
static inline double phip_to_thew(double phi) { return (PI / 2) - phi; }
static inline double lamp_to_phiw(double lambda) {  return lambda; }
static inline double thew_to_phip(double theta) { return (PI / 2) - theta; }
static inline double phiw_to_lamp(double phiw) { return phiw; }

static inline int data_i_ph(const data_iso *data, int i) { return i / data->t; }
static inline int data_i_th(const data_iso *data, int i) { return i % data->t; }

iso_model *newModel(data_iso *data, int lmax);

iso_model *compute_model(int lmax, const char *filename, int npoint);

void modelComputePlm(iso_model *model, const data_iso *data);

void modelComputeCSlm(spherical_model *model, const data_iso *data);

/**========================================================================
 *!                          Coefficients functions
 *========================================================================**/
// Return a newly allocated matrix containing [c_l^m and s_li^m]
// O(lmax^2 * data->N) in time
// the return is a (2 * data->N) x (lmax^2 / 2)
//
// [[c00 c10 c11 ... clm](1)
//  [c00 c10 c11 ... clm](2)
//  [ .      .        . ](.)
//  [ .       .       . ](.)
//  [c00       .     clm](i)
//  [s00 s10 s11 ... clm](1)
//  [s00 s10 s11 ... clm](2)
//  [ .      .        . ](.)
//  [ .       .       . ](.)
//  [ .        .      . ](.)
//  [s00             slm](i)]
//
Matrix_d *compute_mse_coeff(const data_iso *data, const spherical_model *model);

Matrix_d *compute_mse_coeff_clm(const data_iso *data, const spherical_model* model);

Matrix_d *compute_mse_coeff_slm(const data_iso *data, const spherical_model* model);

double compute_pcs_i(const spherical_model *model, const Matrix_d *clm, const Matrix_d *slm, int i);

Matrix_d *compute_pcs_indices(const iso_model *iso, const Matrix_i *indices);

Matrix_d *compute_pcs(const data_iso *data, const spherical_model *model, const Matrix_d *clm, const Matrix_d *slm);

double compute_gradient_clm(const data_iso *data, const Matrix_d *clm, const Matrix_d *pcs, int l, int m);

double compute_gradient_slm(const data_iso *data, const Matrix_d *slm, const Matrix_d *pcs, int l, int m);

double compute_gradient_clm_points(const iso_model *iso, const Matrix_d *clm, const Matrix_d *pcs, int l, int m, const Matrix_i *indices);

double compute_gradient_slm_points(const iso_model *iso, const Matrix_d *slm, const Matrix_d *pcs, int l, int m, const Matrix_i *indices);

double compute_mse(const iso_model *model);

#endif // GEODESY_MODEL_H