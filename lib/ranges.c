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

bool models_equal(const SphericalModel *lhs, const SphericalModel *rhs) {

    if (lhs->lmax != rhs->lmax) return false;
    if (lhs->ll != rhs->ll) return false;

    for (int i = 0; i < lhs->ll; i++) {
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

// Compute the range[a, b] of C and L coefficients
range *modelComputeRange(int a, int b, const data_iso *data, const Precomp *precomp) {

    // if (this_rank == 0)
    range *r = newRange(a, b, NULL, NULL); // Coefficients initialized as zeros(1, ((b + 1)(b + 2) - a(a + 1)) / 2)

    // We want to iterate starting at a and ending at b

    // Matrix_d *C_lm = r->C_lm, *S_lm = r->S_lm;
    // Matrix_d *P_lm_th = precomp->Plm_th;

    printf("[modelComputeRange] Computing coefficient range [%d, %d] in serial\n", a, b);

    Clock *clock = Clock_new();
    Clock_tic(clock);

    cs_pair cs;

    // A single iteration of this loop can be abstracted away as:
    // computeCSPair(l, m, P_lm_th, precomp)

    for (int l = a, i = 0; l <= b; l++) {
        for (int m = 0; m <= l; m++) {
            
            // printf("Computing [%d]th coefficient\n", i);

            cs = computeCSPairAlt(l, m, data, precomp);
            vecset_d(r->C_lm, i, cs.c);
            vecset_d(r->S_lm, i, cs.s);
            i++;
        }
    }

    Clock_toc(clock);
    // double time = elapsed_time(clock);
    free(clock);

    // printf("[modelComputeRange] took %lf s\n", time);

    return r;

}

// Compute the range[a, b] of C and S coefficients
range *modelComputeRangeMPI(int a, int b, const data_iso *data, const Precomp *precomp, int world_size, int this_rank) {

    // range *r = newRange(1, 2, NULL, NULL);
    range *r = newRange(a, b, NULL, NULL);

    // if (this_rank == 0) {
        // freeRange(r);
        // r = newRange(a, b, NULL, NULL); // Coefficients initialized as zeros(1, ((b + 1)(b + 2) - a(a + 1)) / 2)
    // }

    // We want to iterate starting at a and ending at b

    // Matrix_d *C_lm = r->C_lm, *S_lm = r->S_lm;
    // Matrix_d *P_lm_th = precomp->Plm_th;

    if (this_rank == 0) 
        printf("[modelComputeRangeMPI] Computing coefficient range [%d, %d] with %d threads\n", a, b, world_size);

    Clock *clock = Clock_new();
    Clock_tic(clock);

    // cs_pair cs;

    // Total work is the cardinality of my range
    int total_workload = r->card / 2;

    // Go ahead and compute the workload needed by this rank
    int this_workload = compute_workload(total_workload, world_size, this_rank);


    Matrix_i *start_end = compute_startend_array(total_workload, world_size);

    int i_start = vecat_i(start_end, this_rank) + 1;
    int i_end = vecat_i(start_end, this_rank + 1);

    // Allocate the space for Clm and Slm matrices
    Matrix_d *Clm = Matrix_new_d(1, this_workload);
    Matrix_d *Slm = Matrix_new_d(1, this_workload);

    MPI_Barrier(MPI_COMM_WORLD);

    // Here we need to assert the validity of data and precomp
    assert(data != NULL);
    assert(precomp != NULL);
    assert(Clm != NULL);
    assert(Slm != NULL);

    MPI_Barrier(MPI_COMM_WORLD);

    int range_start = r->a * (r->a + 1) / 2;
    // int a_pos = PT(a, 0);

    for (int i = i_start, count = 0; i <= i_end; i++, count++) {
        lm_pair pair = i_to_lm(i + range_start);
        cs_pair cs = computeCSPairAlt(pair.l, pair.m, data, precomp);
        vecset_d(Clm, count, cs.c);
        vecset_d(Slm, count, cs.s);
    }

    // MPI_ORDERED( 
    //     printf("[%d] workload: %d, total_work: %d\n", this_rank, this_workload, total_workload);
    //     printf("[%d] computed %lu coefficients in [%d, %d] with i \\in [%d, %d]\n", this_rank, Matrix_size_d(Clm) * 2, r->a, r->b, i_start, i_end);
    // )

    // Vector_print_head_d(Clm, 20);
    // Vector_print_head_d(Slm, 20);

    Matrix_i *recvcounts = compute_workload_array(total_workload, world_size);
    Matrix_i *displacements = compute_displacements(total_workload, world_size);

    // if (this_rank == 0) {
    //     Matrix_print_i(recvcounts);
    //     Matrix_print_i(displacements);
    // }

    MPI_Barrier(MPI_COMM_WORLD);

    // C coefficients
    MPI_Gatherv(Clm->data, 
        this_workload, 
        MPI_DOUBLE, 
        r->C_lm->data, 
        recvcounts->data, 
        displacements->data, 
        MPI_DOUBLE, 
        0, 
        MPI_COMM_WORLD);

    // S coefficients
    MPI_Gatherv(Slm->data, 
        this_workload, 
        MPI_DOUBLE, 
        r->S_lm->data, 
        recvcounts->data, 
        displacements->data, 
        MPI_DOUBLE, 
        0, 
        MPI_COMM_WORLD);

    // if (this_rank == 0) {
    //     printf("===== Root received Slm->data and Clm->data =====\n");
    // }

    // if (this_rank == 0) {

    //     Matrix_print_d(r->C_lm);
    //     Matrix_print_d(r->S_lm);
    // }

    // Clock_toc(clock);
    // double time = elapsed_time(clock);
    // free(clock);
    
    // Matrix_free_i(recvcounts);
    // Matrix_free_i(displacements);
    // Matrix_free_d(Clm);
    // Matrix_free_d(Slm);
    // Matrix_free_i(start_end);
    // free(clock);

    if (this_rank == 0) {
        return r;
    } else {
        freeRange(r);
        return NULL;
    }

    return r;

}

SphericalModel *SphericalModelAddRange(const SphericalModel *model, const range *r) {

    range *r_model = SphericalModelToRange(model);
    range *r_union = ranges_union(r_model, r);

    // Transfer the matrices from this union to the model.

    SphericalModel *new_model = newSphericalModel(r->b);

    new_model->C_lm = r_union->C_lm;
    new_model->S_lm = r_union->S_lm;

    r_union->C_lm = NULL;
    r_union->S_lm = NULL; // don't manage memory of matrices anymore

    freeRange(r_model);
    freeRange(r_union);

    return new_model;

}

