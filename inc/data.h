#ifndef GEODESY_DATA_H
#define GEODESY_DATA_H

/**========================================================================
 * ?                          data.h
 * @brief   : Read and write functionalities to interoperate with data
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-15
 *========================================================================**/
#include "harmonics.h"
#include "ejovo.h"


// Representation of a sequence of data points following the physics convention
// ISO-80000-2:2019
// th \in [0, 2pi]
// ph \in [0, pi]
typedef struct data_points_iso {
    int npoint;
    double *th;        // polar angle, angle with respect to the polar (z) axis
    double *ph;        // azimuthal angle
    double *r;         // radius, corresponds to f(th, ph)
} data_points_iso;

// condensed version of the data
typedef struct data_iso {
    int N;
    int t;
    int p;
    Matrix_d *th; // only store t values of theta
    Matrix_d *ph; // only store p values of ph
    // Matrix_d *r;  // store N values
    int16_t *r;
    double dt;
    double dp;
    // Matrix_i *l_indices; // Store linear indices of l and m
    // Matrix_i *m_indices;
    const char *size_dataset;
} data_iso;

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

// Free a data pointer that was allocated with malloc
void free_data_iso(data_iso *data);

// Print a summary of the data to the console, printing the first n elements of r 
void head_data(const data_iso *data);

#endif // GEODESY_DATA_H