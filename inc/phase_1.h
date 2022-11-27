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
#include <omp.h>
#include <mpi.h>

#include "constants.h"
#include "data.h"
#include "harmonics.h"
#include "learn.h"
#include "model.h"
#include "mpi_util.h"


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
Matrix_d *read_binary_plm(int lmax, int n_theta, const char *binary_filename, bool __log);

// Return a matrix containing only the P_lm corresponding to l = L
// Return value is a n_theta x (L + 1) matrix
Matrix_d *read_binary_plm_l(int lmax, int l, int n_theta, const char *binary_filename);

// Need a RANGE function that returns a precomputed set of L's between a certain range
Matrix_d *read_binary_plm_range(int a, int b, int n_th, const char *binary_filename, bool __log);

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
SphericalModel *newSphericalModel(int lmax);

double modelComputeCSlmPrecomp(SphericalModel *model, const data_iso *data, const Precomp *precomp);

// Compute all the coefficients Clm and Slm for a model using OMP parallelization at the level of integration. 
// The for loop to compute a single Clm or Slm value is split across omp_get_num_threads() threads. We find that this
// parallelization scheme is incredibly inefficient and actually increases the runtime as we increase the number of 
// threads
double modelComputeCSlmPrecompOMP(SphericalModel *model, const data_iso *data, const Precomp *precomp);

// Compute all the coefficients Clm and Slm for a model using OMP parallelization at the coefficient level.
// The foor loop used to iterate across coefficient computation is split across omp_get_num_threads()
// This implementation is incredibly efficient and leads to a scalability of nearly 100%
double modelComputeCSlmPrecompOMP2(SphericalModel *model, const data_iso *data, const Precomp *precomp);

// Identical to modelComputeCSlmPrecompOMP2 except with the option to specify the number of threads
double modelComputeCSlmPrecompOMP2Threads(SphericalModel *model, const data_iso *data, const Precomp *precomp, int nthreads);

SphericalModel *buildSphericalModelMPI(const data_iso *data, int lmodel, int lbin, const char *coeff_file_bin, bool recompute, bool from, int a, int world_size, int this_rank);

// Write the coefficients of a spherical model to CSV
// Use the prefix to indicate what data set the spherical model was trained on
// pass "" as type if you don't want to specify
void SphericalModelToTXT(const SphericalModel *model, const char *type);

// Write the coefficients of a spherical model to a binary file
// Use the prefix to indicate what data set the spherical model was trained on
void SphericalModelToBIN(const SphericalModel *model, const char *type);

// Load SphericalModel from a binary file
SphericalModel *loadSphericalModel(const char *bin_in, int lmax);

void freeSphericalModel(SphericalModel *model);

typedef struct {

    double c;
    double s;

} cs_pair;

// Compute AND SET the Clm and Slm pair for a given l and m
static inline void computeCSPair(int l, int m, const data_iso *data, const Precomp *precomp, Matrix_d *C_lm, Matrix_d *S_lm) {

    // Now that I have the Associated legendre functions and the data efficiently loaded, let's write
    // the code to approximate the integrals
    double c_integral = 0;
    double s_integral = 0;
    int count = 0;

    const Matrix_d *P_lm_th = precomp->Plm_th;

    for (int i = 0; i < data->N; i++) {

        int i_th = i % data->t;
        int i_ph = i / data->t;

        // use the midpoint formula
        c_integral += data->r[i] * matat_d(P_lm_th, i_th, PT(l, m)) * matat_d(precomp->cosmph, m, i_ph) * vecat_d(precomp->sinth, i_th);
        s_integral += data->r[i] * matat_d(P_lm_th, i_th, PT(l, m)) * matat_d(precomp->sinmph, m, i_ph) * vecat_d(precomp->sinth, i_th);

        count ++;
    }

    c_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);
    s_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);

    //! WARNING this is BLACK MAGIC. I don't know MATHEMATICALLY why this even works... but it does!
    if (m != 0) {
        c_integral *= 4;
        s_integral *= 4;
    }

    vecset_d(C_lm, PT(l, m), c_integral);
    vecset_d(S_lm, PT(l, m), s_integral);

}

// Compute the Clm and Slm coefficients and return the value (instead of setting at a guess of the position) 
static inline cs_pair computeCSPairAlt(int l, int m, const data_iso *data, const Precomp *precomp) {

    // Now that I have the Associated legendre functions and the data efficiently loaded, let's write
    // the code to approximate the integrals
    double c_integral = 0;
    double s_integral = 0;
    int count = 0;

    const Matrix_d *P_lm_th = precomp->Plm_th;

    for (int i = 0; i < data->N; i++) {

        int i_th = i % data->t;
        int i_ph = i / data->t;

        // use the midpoint formula
        c_integral += data->r[i] * matat_d(P_lm_th, i_th, PT(l, m)) * matat_d(precomp->cosmph, m, i_ph) * vecat_d(precomp->sinth, i_th);
        s_integral += data->r[i] * matat_d(P_lm_th, i_th, PT(l, m)) * matat_d(precomp->sinmph, m, i_ph) * vecat_d(precomp->sinth, i_th);

        count ++;
    }

    c_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);
    s_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);

    //! WARNING this is BLACK MAGIC. I don't know MATHEMATICALLY why this even works... but it does!
    if (m != 0) {
        c_integral *= 4;
        s_integral *= 4;
    }

    cs_pair cs = { .c = c_integral, .s = s_integral };
    return cs;
}

// Compute AND SET the Clm and Slm pair for a given l and m
// Change how we fill the C_lm and S_lm matrices.
static inline void computeCSPairMPI(int l, int m, const data_iso *data, const Precomp *precomp, Matrix_d *C_lm, Matrix_d *S_lm, int i_start) {

    // Now that I have the Associated legendre functions and the data efficiently loaded, let's write
    // the code to approximate the integrals
    double c_integral = 0;
    double s_integral = 0;
    int count = 0;

    const Matrix_d *P_lm_th = precomp->Plm_th;

    for (int i = 0; i < data->N; i++) {

        int i_th = i % data->t;
        int i_ph = i / data->t;

        // use the midpoint formula
        c_integral += data->r[i] * matat_d(P_lm_th, i_th, PT(l, m)) * matat_d(precomp->cosmph, m, i_ph) * vecat_d(precomp->sinth, i_th);
        s_integral += data->r[i] * matat_d(P_lm_th, i_th, PT(l, m)) * matat_d(precomp->sinmph, m, i_ph) * vecat_d(precomp->sinth, i_th);

        count ++;
    }

    c_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);
    s_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);

    //! WARNING this is BLACK MAGIC. I don't know MATHEMATICALLY why this even works... but it does!
    if (m != 0) {
        c_integral *= 4;
        s_integral *= 4;
    }

    vecset_d(C_lm, PT(l, m) - i_start, c_integral);
    vecset_d(S_lm, PT(l, m) - i_start, s_integral);

}

// Compute AND SET the Clm and Slm pair for a given l and m using OpenMP
static inline void computeCSPairOMP(int l, int m, const data_iso *data, const Precomp *precomp, Matrix_d *C_lm, Matrix_d *S_lm, Matrix_d *P_lm_th) {

    // Now that I have the Associated legendre functions and the data efficiently loaded, let's write
    // the code to approximate the integrals
    double c_integral = 0, c_partial = 0;
    double s_integral = 0, s_partial = 0;
    int count = 0;

    #pragma omp parallel private(c_partial, s_partial) shared(c_integral, s_integral)
    {

        #pragma omp for 
        for (int i = 0; i < data->N; i++) {

            int i_th = i % data->t;
            int i_ph = i / data->t;

            // use the midpoint formula
            c_partial += data->r[i] * matat_d(P_lm_th, i_th, PT(l, m)) * matat_d(precomp->cosmph, m, i_ph) * vecat_d(precomp->sinth, i_th);
            s_partial += data->r[i] * matat_d(P_lm_th, i_th, PT(l, m)) * matat_d(precomp->sinmph, m, i_ph) * vecat_d(precomp->sinth, i_th);

            count ++;
        }

        #pragma omp critical
        {
            c_integral += c_partial;
            s_integral += s_partial;
        }
    }

    c_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);
    s_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);

    //! WARNING this is BLACK MAGIC. I don't know MATHEMATICALLY why this even works... but it does!
    if (m != 0) {
        c_integral *= 4;
        s_integral *= 4;
    }

    vecset_d(C_lm, PT(l, m), c_integral);
    vecset_d(S_lm, PT(l, m), s_integral);

}

/**========================================================================
 *!                           MPI Implementations
 *========================================================================**/
// Top-level function to compute all of the coefficients using MPI parallelization techniques
// Only the root process will end up filling the model.
double modelComputeCSlmPrecompMPI(SphericalModel *model, const data_iso *data, const Precomp *precomp, int world_size, int this_rank);

// Allocate space for and load (or compute) a spherical model of degree lmodel
SphericalModel *buildSphericalModel(const data_iso *data, int lmodel, const char *coeff_file_bin, bool recompute, bool from, int a);

double modelComputeCSlmPrecompAlt(SphericalModel *model, const data_iso *data, const Precomp *precomp);

#endif //GEODESY_PHASE_1_H