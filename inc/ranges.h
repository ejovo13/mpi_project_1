#ifndef GEODESY_RANGES_H
#define GEODESY_RANGES_H
/**========================================================================
 * ?                          ranges.h
 * @brief   : Implementation of ranges as discussed in the project report
 * @details : C_[a, b] = {C_a^0, ..., C_a^a, ..., C_b^0, ..., C_b^b} and is
 * stored in a laplace_range_t
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-11-23
 *========================================================================**/
#include "ejovo.h"
#include "model.h"
#include "data.h"
#include "phase_1.h"
#include "phase_2.h"
#include "phase_3.h"

typedef struct laplace_range_t {

    int a;
    int b;
    Matrix_d* C_lm;
    Matrix_d* S_lm;
    int card; // Total number of coefficients, card := |C_lm| + |S_lm| = (b + 1)(b + 2) - (a)(a + 1)

} range;

static inline void freeRange(range *r) {

    if (!r) return;

    if (r->C_lm) 
        Matrix_free_d(r->C_lm);

    if (r->S_lm)
        Matrix_free_d(r->S_lm);

    free(r);
}

// Allocate space for and return a new range, cloning the input matrices.
// If C_lm or S_lm are null, they will be replaced with a zeros vector
range* newRange(int a, int b, const Matrix_d* C_lm, const Matrix_d* S_lm);

// Contribute t
// void addToRange

// void add_to

range * ranges_union(const range *lhs, const range *rhs);

// Check if the ranges store the IDENTICAL information in their C_lm and S_lm vectors
bool ranges_equal(const range *lhs, const range *rhs);

// Extract the range[a, b] from a computed model
range *SphericalModelExtractRange(const SphericalModel *model, int l_0, int l_f);

// Check to see if the intervals of lhs and rhs are exactly equal
bool same_intervals(const range *lhs, const range *rhs);

void printRange(const range *r);

range *SphericalModelToRange(const SphericalModel *model);

range *as_range(const SphericalModel *model);

void printRangeCoeff(const range *r);

range *modelComputeRange(int a, int b, const data_iso *data, const Precomp *precomp);

range *modelComputeRangeMPI(int a, int b, const data_iso *data, const Precomp *precomp, int world_size, int this_rank);

// Add a range to a spherical model, returning a new model
SphericalModel *SphericalModelAddRange(const SphericalModel *model, const range *r);

bool models_equal(const SphericalModel *lhs, const SphericalModel *rhs);

#endif // GEODESY_RANGES_H