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

#include "ejovo.h"

/* representation of a sequence of data points */
struct data_points {
	int npoint;
	double *phi;
	double *lambda;
	double *V;
};

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
    Matrix_d *r;  // store N values
    double dt;
    double dp;
    Matrix_i *l_indices; // Store linear indices of l and m
    Matrix_i *m_indices;
} data_iso;



#endif