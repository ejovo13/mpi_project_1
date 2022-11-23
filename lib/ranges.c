#include "ranges.h"

// Allocate space for and return a new range, cloning the input matrices.
// If C_lm or S_lm are null, they will be replaced with a zeros vector
range* newRange(int a, int b, const Matrix_d* C_lm, const Matrix_d* S_lm) {

    range* r = (range *) malloc(sizeof(*r));
    r->a = a;
    r->b = b;
    r->card = (b + 1) * (b + 2) - a * (a + 1);

    if (C_lm == NULL) {
        r->C_lm = Matrix_new_d(1, r->card / 2);
    } else {
        r->C_lm = Matrix_clone_d(C_lm);
    }
    if (S_lm == NULL) {
        r->S_lm = Matrix_new_d(1, r->card / 2);
    } else {
        r->S_lm = Matrix_clone_d(S_lm);
    }

    return r;
}

// Contribute t
// void addToRange

// void add_to

range * ranges_union(const range *lhs, const range *rhs) {
    // We don't check overlap between lhs.b and rhs.a (maybe we should) 
    // compute the new bounds, so that we know how big to make the final matrices
    assert(lhs->a < rhs->a);
    // assert(rhs->a > lhs->b);  // this is not necessary

    int new_a = lhs->a < rhs->a ? lhs->a : rhs->a;
    int new_b = lhs->b > rhs->b ? lhs->b : rhs->b;

    // Initialize all the values equal to 0, then add them in later
    range *out = newRange(new_a, new_b, NULL, NULL);

    // We don't iterate along the elements of the output, we iterate along
    // the two input ranges

    // Compute the index in the final range (new_a, new_b)
    // the value of PT(new_a, 0).  
    int last_index_lhs = PT(lhs->a, 0);
    
    // assert(start_index == last_index_lhs);
    for (int i = 0; i < lhs->card / 2; i++) {
        *vecacc_d(out->C_lm, i) = vecat_d(lhs->C_lm, i);
        *vecacc_d(out->S_lm, i) = vecat_d(lhs->S_lm, i);
    }

    // For two overlapping intervals, we want to 

    int last_index_rhs = PT(rhs->a, 0);
    for (int i = 0; i < rhs->card / 2; i++) {
        *vecacc_d(out->C_lm, last_index_rhs - last_index_lhs + i) = vecat_d(rhs->C_lm, i);
        *vecacc_d(out->S_lm, last_index_rhs  - last_index_lhs+ i) = vecat_d(rhs->S_lm, i);
    }

    return out;

}

// Check if the ranges store the IDENTICAL information in their C_lm and S_lm vectors
bool ranges_equal(const range *lhs, const range *rhs) {
    if (!same_intervals(lhs, rhs)) return false;

    for (int i = 0; i < lhs->card / 2; i++) {
        if (lhs->C_lm->data[i] != rhs->C_lm->data[i]) return false;
        if (lhs->S_lm->data[i] != rhs->S_lm->data[i]) return false;
    }

    return true;
}

// Extract the range[a, b] from a computed model
range *SphericalModelExtractRange(const SphericalModel *model, int l_0, int l_f) {

    printf("trying to extract range [%d, %d]\n", l_0, l_f);

    range *out = newRange(l_0, l_f, NULL, NULL);

    int init_pos = ((l_0) * (l_0 + 1) / 2);

    printf("Out->card %d\n", out->card);

    printf("init_pos: %d\n", init_pos);

    for (int i = 0; i < out->card / 2; i++) {
        *vecacc_d(out->C_lm, i) = vecat_d(model->C_lm, init_pos + i);
        *vecacc_d(out->S_lm, i) = vecat_d(model->S_lm, init_pos + i);
    }

    return out;
}

// Check to see if the intervals of lhs and rhs are exactly equal
bool same_intervals(const range *lhs, const range *rhs) {
    return (lhs->a == rhs->a && lhs->b == rhs->b);
}

void printRange(const range *r) {
    printf("range[%d, %d], .card = %d\n", r->a, r->b, r->card);
}

void printRangeCoeff(const range *r) {
    printf("range[%d, %d], .card = %d\n", r->a, r->b, r->card);
    Matrix_print_d(r->C_lm);
    Matrix_print_d(r->S_lm);
}

range *SphericalModelToRange(const SphericalModel *model) {
    return newRange(0, model->lmax, model->C_lm, model->S_lm);
}

range *as_range(const SphericalModel *model) {
    return SphericalModelToRange(model);
}

