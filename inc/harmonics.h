#ifndef GEODESY_H
#define GEODESY_H
#include <assert.h>
#include "ejovo.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
#endif

#ifndef TWO_PI
#define TWO_PI (2.0 * PI)
#endif

#ifndef HALF_PI
#define HALF_PI (PI / 2.0)
#endif

#ifndef SPHERICAL_MODEL_CONSTANTS
#define SPHERICAL_MODEL_CONSTANTS

    // ranges of th and ph
    #define TH_0 0
    #define TH_F PI
    #define PH_0 0
    #define PH_F TWO_PI

#endif

/* representation of spherical harmonics coefficients */
struct spherical_harmonics {
	int lmax;
	double *CS;
	double *A;
	double *B;
};

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

// extract the value of Clm from a model



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


/* these 3 functions help compute indices in arrays containing
   triangular matrices */
static inline int PT(int l, int m)   /* classic triangular w/ diagonal */
{
	return m + l * (l + 1) / 2;
}

static inline int CT(int l, int m)
{
	assert(m <= l);
	return m + l * l;
}

static inline int ST(int l, int m)   /* all Sx0 are missing */
{
	assert(m <= l);
	assert(1 <= m);
	return m + (l+1) * l;
}

static inline int PLM(int l, int m, int th, int nrows) {
    int index = PT(l, m);
    return th * nrows + index; 
}

// Assume that the model has already been initialized
static inline double model_C(const spherical_model *model, int l, int m) {
    return model->C_lm->data[PT(l, m)];
}

static inline void model_set_C(spherical_model *model, int l, int m, double val) {
    model->C_lm->data[PT(l, m)] = val;
}

static inline void model_C_plus(spherical_model *model, int l, int m, double addition) {
    model->C_lm->data[PT(l, m)] += addition;
}

static inline double model_S(const spherical_model *model, int l, int m) {
    return model->S_lm->data[PT(l, m)];
}

// conversion functions from iso to project
// project to wolfram
static inline double phip_to_thew(double phi) {
    return (PI / 2) - phi;
}

static inline double lamp_to_phiw(double lambda) {
    // return lambda + PI;
    return lambda;
}

static inline double thew_to_phip(double theta) {
    return (PI / 2) - theta;
}

static inline double phiw_to_lamp(double phiw) {
    // return phiw - PI;
    return phiw;
}

static inline int data_i_ph(const data_iso *data, int i) {
    return i / data->t;
}

static inline int data_i_th(const data_iso *data, int i) {
    return i % data->t;
}


/**========================================================================
 *!                           Opeartions on Models
 *========================================================================**/
// top-level functions
double compute_mse(const iso_model *model);


// What do I need to create a new model? A set of data and an lmax
void modelComputePlm(iso_model *model, const data_iso *data);
void modelComputeCSlm(spherical_model *model, const data_iso *data);
// spherical_model *newModel(const data_iso *data, int lmax);
iso_model *newModel(data_iso *data, int lmax);
Matrix_d *compute_gradient(const data_iso *data, const spherical_model *model);
void adjust_parameters(const Matrix_d *grad, spherical_model *model, double alpha);
double compute_average_error(const iso_model *iso);
void writeModel(const spherical_model *model, const data_iso *data, const char *prefix);

double modelPredict(const iso_model *iso, double theta, double phi);
double modelPredictDataPoint(const iso_model *iso, int i);
Matrix_d *modelPredictDataPoints(const iso_model *iso, Matrix_i *indices);
Matrix_d *modelPredictData(const iso_model *iso);
Matrix_d *modelPredictN(const iso_model *iso, int __n);

double estimate_mse(const iso_model *iso, int n);
// Retrieve the value of P_l^m(cos \theta_i)
// 
// P_lm_th is a n_theta x ll matrix 
//                                                                     theta_i
static inline double model_P(spherical_model *model, int l, int m, int thi) {
    return matat_d(model->P_lm_th, thi, PT(l, m));
} 


/* wall-clock seconds (with high precision) elapsed since some given point in the past */
double wtime();

/* represent n in <= 8 char, in a human-readable way  */
void human_format(char * target, long n);

/* allocates memory for self */
void setup_spherical_harmonics(int lmax, struct spherical_harmonics *self);

/*
 * Load spherical harmonics from a file.  File format:
 * each line contains 2 integers and 2 floating-point numbers separated by tabs.
 * Each line contains (l, m, c, s) with 0 <= m <= l <= lmax.
 * If m == 0, then s == 0.
 */
void load_spherical_harmonics(const char *filename, int lmax, struct spherical_harmonics *self);

/*
 * Load data points from a file into "self".  File format:
 * each line contains 3 floating-point numbers separated by tabs.
 * Each line contains (lambda, phi, V) where V == f(phi, lambda).
 */
void load_data_points(const char *filename, int npoint, struct data_points *self);


data_points_iso *load_data_points_iso(const char *filename, int npoint);

data_iso *load_data_iso(const char *filename, int t, int p);

void write_iso(const data_iso *data, const char *filename);

void print_npoints(const data_points_iso* data, int n);

void write_npoints(const data_points_iso* data, int npoints, const char *filename);

/*
 * Compute all the (fully normalized) Associated Legendre function of degree <= lmax.
 * On exit, P_{l,m} (with 0 <= m <= l <= lmax) can be found in P[PT(l, m)].
 * P must be preallocated of size (lmax + 1) * (lmax + 2) / 2.
 */
void computeP(const struct spherical_harmonics *self, double *P, double sinphi);

/* phi   : elevation angle (latitude) in degrees north of the equator plane, 
           in the range -PI/2 <= phi <= PI/2 
   lambda: azimuth angle (longitude) measured in degrees east or west from some
           conventional reference meridian.
   P must be previously evaluated with a call to computeP(self, P, sin(phi));
*/
double evaluate(const struct spherical_harmonics *self, const double *P, double lambda);

/**========================================================================
 *!                           Gradient and coefficients functions
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

Matrix_d *compute_pcs(const data_iso *data, const spherical_model *model, const Matrix_d *clm, const Matrix_d *slm);

double compute_gradient_clm(const data_iso *data, const Matrix_d *clm, const Matrix_d *pcs, int l, int m);

double compute_gradient_slm(const data_iso *data, const Matrix_d *slm, const Matrix_d *pcs, int l, int m);

Matrix_d *modelPredictData(const iso_model *iso);

double modelPredictDataPoint(const iso_model *iso, int i);

// Compute gradient in the form 
// [dMSE/dc00 dMSE/dc10 ... dMSE/dclm
//  dMSE/ds00 dMSE/ds10 ... dMSE/dslm]
// Assume model has C_lm and S_lm and P_lm_th
Matrix_d *compute_gradient(const data_iso *data, const spherical_model *model);

Matrix_d *compute_stochastic_gradient(const iso_model *iso, int n);


// Assume that a model object has been fully initialized
double compute_mse(const iso_model *iso);

double compute_average_error(const iso_model *iso);

// Given a gradient matrix and a model, tweak the values of 
// C_lm to improve the score
// gradient has size (N * 2) x ((lmax + 1)(lmax + 2)/2)
void adjust_parameters(const Matrix_d *grad, spherical_model *model, double alpha);

#endif