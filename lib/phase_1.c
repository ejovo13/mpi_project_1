#include "geodesy.h"

void write_binary_plm(int lmax, const Matrix_d *th, const char *binary_file_out) {

    const int LL = (lmax + 1) * (lmax + 2) / 2;

    // This is hardcoded and really a bad idea.
    // should adjust the function to take in the full name of the input 
    // file
    // char filename[100] = {0};
    // sprintf(filename, "ETOPO1_%s_P%d.bin", dataset, lmax);

    // Allocate space for a single row of the matrix

    // initialize coefficients a and b
    struct spherical_harmonics sph_model;
    setup_spherical_harmonics(lmax, &sph_model);

    // Open up the binary file to write
    FILE *bin = fopen(binary_file_out, "wb");

    int n = Matrix_size_d(th);
    int n_bytes_per_row = sizeof(double) * LL;

    Clock *clock = Clock_new();


        Clock_tic(clock);


    Matrix_d *plm = Matrix_new_d(1, LL);

    for (int i = 0; i < n; i++) {
        // Now compute the p values
        computeP(&sph_model, plm->data, cos(th->data[i]));
        fwrite(plm->data, n_bytes_per_row, 1, bin);
    }


        Clock_toc(clock);

    fclose(bin);
    printf("[write_binary_plm] Wrote to binary file: %s in %lfs\n", binary_file_out, elapsed_time(clock));

    free(clock);
}

Matrix_d *read_binary_plm(int lmax, int n_th, const char *binary_filename, bool __log) {

    const size_t LL = (lmax + 1) * (lmax + 2) / 2;

    if (__log) 
        printf("[read_binary_plm] LL: %lu, n_theta: %d\n", LL, n_th);

    // Allocate a Matrix to store the contents
    Matrix_d *P_lm_th = Matrix_new_d(n_th, LL);

    if (__log)
        printf("[read_binary_plm] Allocated new matrix of size %lu x %lu\n", P_lm_th->nrows, P_lm_th->ncols);

    size_t n_bytes_per_row = sizeof(double) * LL;

    FILE *bin = fopen(binary_filename, "rb");

    // now fill up this matrix, row by row
    for (int i = 0; i < n_th; i++) {
        fread(matacc_d(P_lm_th, i, 0), n_bytes_per_row, 1, bin);
    }

    fclose(bin);

    if (__log)
        printf("[read_binary_plm] Read from binary file: %s\n", binary_filename);

    return P_lm_th;

}

// Load only plm from l = a to b
Matrix_d *read_binary_plm_range(int a, int b, int n_th, const char *binary_filename, bool __log) {

    const size_t card = ((b + 1) * (b + 2) - (a) * (a + 1))/ 2;

    if (__log) 
        printf("[read_binary_plm_range] [a, b]: [%d, %d], card: %lu, n_theta: %d\n", a, b, card, n_th);

    // Allocate a Matrix to store the contents
    Matrix_d *P_lm_th = Matrix_new_d(n_th, card);

    if (__log)
        printf("[read_binary_plm] Allocated new matrix of size %lu x %lu\n", P_lm_th->nrows, P_lm_th->ncols);

    size_t n_bytes_per_row = sizeof(double) * card;
    size_t n_bytes_before_a = sizeof(double) * (a * (a + 1) / 2);

    FILE *bin = fopen(binary_filename, "rb");

    if (!bin) {
        errx(2, "[read_binary_plm_range] failed to open %s; exiting\n", binary_filename);
    }

    // now fill up this matrix, row by row
    for (int i = 0; i < n_th; i++) {
        fseek(bin, n_bytes_before_a, SEEK_CUR); // skip the first (a - 1) computed values for each row
        fread(matacc_d(P_lm_th, i, 0), n_bytes_per_row, 1, bin);
    }

    fclose(bin);

    if (__log)
        printf("[read_binary_plm_range] Read from binary file: %s\n", binary_filename);

    return P_lm_th;

}

Matrix_d *read_binary_plm_l(int lmax, int l, int n_theta, const char *binary_filename) {

    const size_t LL = (lmax + 1) * (lmax + 2) / 2; // total number of elements per row in binary file
    const size_t l_middle = (l + 1) * (l + 2) / 2; // useful number to comput the trailing number of elements
    const size_t l_preceding = (l) * (l + 1) / 2;
    const size_t l_trailing = LL - l_middle;

    // We will retrieve L columns from the binary file.
    // First thing that we need to do is compute an offset of bytes
    // that we will seek to every time we try to read

    // For example for L = 0, the offset should be 0
    // For L = 3, we need to pass P00 P10 P11 P20 P21 P22
    const size_t row_offset_bytes = l_preceding * sizeof(double);
    const size_t bytes_per_row = (l + 1) * sizeof(double);
    const size_t row_trailing_bytes = l_trailing * sizeof(double);

    // allocate proper amount of space for the matrix, which will be n_theta x L
    Matrix_d *P_L_th = Matrix_new_d(n_theta, l + 1);

    FILE *bin = fopen(binary_filename, "rb");

    Clock *clock = Clock_new();

    Clock_tic(clock);
        // Now fill the matrix rowwise
        for (int i = 0; i < n_theta; i++) {
            fseek(bin, row_offset_bytes, SEEK_CUR); // seek to start of the row
            fread(matacc_d(P_L_th, i, 0), bytes_per_row, 1, bin);
            fseek(bin, row_trailing_bytes, SEEK_CUR);        
        }
    Clock_toc(clock);

    fclose(bin);
    printf("Read in %d x %d P_L_th matrix from %s\n", n_theta, l + 1, binary_filename);

    free(clock);

    return P_L_th;
}

// TODO THIS FUNCTION COMPLETELY IGNORE LF. NOT GOOD!!!!
Precomp *newPrecomp(int L0, int LF, int Lmax, const data_iso *data, const char *plm_bin) {

    // Allocate space
    Precomp *precomp = (Precomp *) malloc(sizeof(*precomp));
    if (precomp == NULL) 
        errx(1, "[newPrecomp] Problem allocating Precomp\n");

    precomp->L0 = L0;
    precomp->LF = LF;
    
    // Check if file exists. If it doesn't, precompute the values
    FILE *bin = fopen(plm_bin, "rb");
    if (bin == NULL) {
        // then the file doesn't exist, go ahead and create it.
        write_binary_plm(Lmax, data->th, plm_bin);
    }

    // Now retrieve the proper values
    // Ideally I should only be loading L0 - LF into memory
    precomp->Plm_th = read_binary_plm(Lmax, data->t, plm_bin, false);

    // Matrix_print_d(precomp->Plm_th);

    precomp->sinth  = Matrix_new_d(1, data->t);
    // m goes from 0 to L
    precomp->cosmph = Matrix_new_d(Lmax + 1, data->p);
    precomp->sinmph = Matrix_new_d(Lmax + 1, data->p);


    // And compute sin(\theta), sin(m\theta), and cos(m\theta)
    for (int ith = 0; ith < data->t; ith++) {
        precomp->sinth->data[ith] = sin(data->th->data[ith]);
    }

    for (int m = 0; m <= Lmax; m++) {
        for (int iph = 0; iph < data->p; iph++) {
            *matacc_d(precomp->cosmph, m, iph) = cos(m * data->ph->data[iph]);
            *matacc_d(precomp->sinmph, m, iph) = sin(m * data->ph->data[iph]);
        }
    }

    return precomp;
}

void freePrecomp(Precomp *precomp) {

    if (precomp == NULL) return;

    if (precomp->Plm_th) Matrix_free_d(precomp->Plm_th);
    if (precomp->sinth) Matrix_free_d(precomp->sinth);
    if (precomp->sinmph) Matrix_free_d(precomp->sinmph);
    if (precomp->sinmph) Matrix_free_d(precomp->cosmph);

    free(precomp);

}

/**========================================================================
 *!                           Spherical Model
 *========================================================================**/

SphericalModel *newSphericalModel(int lmax) {

    SphericalModel *model = (SphericalModel *) malloc(sizeof(*model));

    if (model == NULL) 
        errx(1, "Problems allocating new SphericalModel\n");

    if (lmax < 0)
        errx(1, "lmax [%d] must be a positive value\n", lmax);

    model->lmax = lmax;
    model->ll = (lmax + 1) * (lmax + 2) / 2;

    // Need to allocate space for my Clm and Slm!!

    model->C_lm = Matrix_new_d(1, model->ll);
    model->S_lm = Matrix_new_d(1, model->ll);

    // Wait a second...
    // // Now let's compute the coefficients (all the coefficients for now)
    // // We are not at the point where we can compute the coefficients individually,
    // // although we _will_ get there

    return model;
}

// Allocate the space for and find the coefficients of a spherical model M of degree l_max
SphericalModel *buildSphericalModel(const data_iso *data, int lmodel, const char *coeff_file_bin, bool recompute, bool from, int a) {

    char plm_bin[100] = {0};
    // Move precomp processing here to remove a parameter in the interface.
    const int lbin = lmodel;
    sprintf(plm_bin, "ETOPO1_%s_P%d.bin", data->size_dataset, lbin);
    Precomp *precomp = newPrecomp(0, lmodel, lbin, data, plm_bin);

    SphericalModel *model = NULL; // Declare our model pointer

    char from_a_bin[100] = {0};
    Clock *clock = Clock_new();

    // Let's first check to see if the binary file already exists. If it exists, load from it.
    FILE *test_open = fopen(coeff_file_bin, "rb");
    if (test_open == NULL || recompute) {

        // double time = modelComputeCSlmPrecomp(model, data, precomp);
        if (recompute) {
            printf("[buildSphericalModel] recomputing coefficients\n");
        } else {
            printf("[buildSphericalModel] %s not found, computing Clm and Slm coefficients\n", coeff_file_bin);
        } 
        double time = 0;

        range *r;

        if (from) {

            sprintf(from_a_bin, "sph_%s_%d.bin", data->size_dataset, a);
            // build a first model up to the a degree
            printf("[buildSphericalModel] First building model from a: %d\n", a);
            model = buildSphericalModel(data, a, from_a_bin, false, false, 0);
            printf("[buildSphericalModel] Now computing range [%d, %d]\n", a, lmodel);
            Clock_tic(clock);
            r = modelComputeRange(a, lmodel, data, precomp);
            Clock_toc(clock);
            model = SphericalModelAddRange(model, r);
            time = elapsed_time(clock);

            // printRange(r);
        } else {
            // Recompute all of the coefficients
            printf("[buildSphericalModel] Computing full set of coefficients\n");
            model = newSphericalModel(lmodel);
            Clock_tic(clock);
            // time = modelComputeCSlmPrecompAlt(model, data, precomp);
            time = modelComputeCSlmPrecomp(model, data, precomp);
            Clock_toc(clock);
            
        }

        printf("[buildSpherialModel] Computed coefficients for L = %d in %lfs\n", lmodel, time);
        SphericalModelToBIN(model, data->size_dataset);

    } else {

        // Model file coeff_file_bin already exists, so load from it
        printf("[buildSphericalModel] Coefficients already present in %s, loading them now\n", coeff_file_bin);
        fclose(test_open);
        model = loadSphericalModel(coeff_file_bin, lmodel);

    }

    freePrecomp(precomp);
    free(clock);    

    return model;

}

// Allocate the space for and find the coefficients of a spherical model M of degree l_max
// I want to sent a copy of the coefficients to EVERY process
SphericalModel *buildSphericalModelMPI(const data_iso *data, int lmodel, int lbin, const char *coeff_file_bin, bool recompute, bool from, int a, int world_size, int this_rank) {

    char plm_bin[100] = {0};
    sprintf(plm_bin, "ETOPO1_%s_P%d.bin", data->size_dataset, lbin);

    // Precomp *precomp

    Precomp *precomp;

    MPI_ORDERED(
        precomp = newPrecomp(0, lmodel, lbin, data, plm_bin);
    )

    SphericalModel *model = NULL; // Declare our model pointer

    char from_a_bin[100] = {0};
    Clock *clock = Clock_new();

    // Let's first check to see if the binary file already exists. If it exists, load from it.
    FILE *test_open = fopen(coeff_file_bin, "rb");
    if (test_open == NULL || recompute) {

        // double time = modelComputeCSlmPrecomp(model, data, precomp);
        MPI_ONCE (
            if (recompute) {
                printf("[buildSphericalModel] recomputing coefficients\n");
            } else {
                printf("[buildSphericalModel] %s not found, computing Clm and Slm coefficients\n", coeff_file_bin);
            } 
        )
        double time = 0;

        range *r;

        // This section computes the model and stores it in model for the root rank
        if (from) {

            sprintf(from_a_bin, "sph_%s_%d.bin", data->size_dataset, a);
            // build a first model up to the a degree
            MPI_ONCE (
                printf("[buildSphericalModel] First building model from a: %d\n", a);
            )

            MPI_Barrier(MPI_COMM_WORLD);
            model = buildSphericalModelMPI(data, a, lbin, from_a_bin, false, false, 0, world_size, this_rank);
            MPI_ONCE (
                printf("[buildSphericalModel] Now computing range [%d, %d]\n", a, lmodel);
            )
            
            Clock_tic(clock);

            // r is only legit in the ROOT PROCESS
            r = modelComputeRangeMPI(a, lmodel, data, precomp, world_size, this_rank);
            Clock_toc(clock);

            if (this_rank == 0) {
                assert(r);
                model = SphericalModelAddRange(model, r);
            }

            time = elapsed_time(clock);

            // // Now save the model so that the other processes can access it.
            // if (this_rank == 0) {
            //     SphericalModelToBIN(model, coeff_file_bin);
            // }

            // MPI_Barrier(MPI_COMM_WORLD);

            // if (this_rank != 0) {
            //     model = loadSphericalModel(coeff_file_bin, lmodel);
            // }

            if (this_rank == 0) {
                freeRange(r);
            }
            // printRange(r);
        } else {
            // Recompute all of the coefficients
            MPI_ONCE(
                printf("[buildSphericalModel] Computing full set of coefficients\n");
            )
            model = newSphericalModel(lmodel);
            Clock_tic(clock);
            // time = modelComputeCSlmPrecompAlt(model, data, precomp);
            time = modelComputeCSlmPrecompMPI(model, data, precomp, world_size, this_rank);
            Clock_toc(clock);
            
        }

        MPI_ONCE (
            printf("[buildSpherialModel] Computed coefficients for L = %d in %lfs\n", lmodel, time);
            SphericalModelToBIN(model, data->size_dataset);
        )

        // There should now be a binary file storing the coefficients. Let's lock up, and then have everybody read from this 
        // file to load their own model. This creates a copy of model for every process.
        MPI_Barrier(MPI_COMM_WORLD);

        if (this_rank != 0) {
            model = loadSphericalModel(coeff_file_bin, lmodel);
        }

        if (this_rank == 0) printf("All models loaded\n");

    } else {

        // Model file coeff_file_bin already exists, so load from it
        MPI_ONCE (
            printf("[buildSphericalModel] Loading coefficients for L = %d\n", lmodel);
        )
        fclose(test_open);
        model = loadSphericalModel(coeff_file_bin, lmodel);

    }

    free(clock);
    freePrecomp(precomp);

    return model;

}

// call computeCSPairAlt to return a { .c, .l } structure instead of setting in place.
double modelComputeCSlmPrecompAlt(SphericalModel *model, const data_iso *data, const Precomp *precomp) {

    const int lmax = model->lmax;
    // const int ll = model->ll;

    printf("[modelComputeCSlmPrecompAlt] Computing coefficients in serial\n");

    Clock *clock = Clock_new();
    Clock_tic(clock);

    // A single iteration of this loop can be abstracted away as:
    // computeCSPair(l, m, P_lm_th, precomp)
    for (int l = 0, i = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++, i++) {
            cs_pair cs = computeCSPairAlt(l, m, data, precomp);

            model->C_lm->data[i] = cs.c;
            model->S_lm->data[i] = cs.s;
        }
    }

    Clock_toc(clock);
    double time = elapsed_time(clock);
    free(clock);

    return time;

}
// Compute the value of Clm and Slm coefficients using the new Precomp data structure
double modelComputeCSlmPrecomp(SphericalModel *model, const data_iso *data, const Precomp *precomp) {

    const int lmax = model->lmax;
    // const int ll = model->ll;

    Matrix_d *C_lm = model->C_lm, *S_lm = model->S_lm;
    // Matrix_d *P_lm_th = precomp->Plm_th;

    printf("[modelComputeCSlmPrecomp] Computing coefficients in serial\n");

    Clock *clock = Clock_new();
    Clock_tic(clock);

    // A single iteration of this loop can be abstracted away as:
    // computeCSPair(l, m, P_lm_th, precomp)
    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            // printf("[modelComputeCSlmPrecomp] computing C[%d,%d]\n", l, m);
            computeCSPair(l, m, data, precomp, C_lm, S_lm);
        }
    }

    Clock_toc(clock);
    double time = elapsed_time(clock);
    free(clock);

    return time;

}

// Compute the value of Clm and Slm coefficients using the new Precomp data structure and
// using a parallel loop with OpenMP
double modelComputeCSlmPrecompOMP(SphericalModel *model, const data_iso *data, const Precomp *precomp) {

    const int lmax = model->lmax;
    // const int ll = model->ll;

    Matrix_d *C_lm = model->C_lm, *S_lm = model->S_lm;
    Matrix_d *P_lm_th = precomp->Plm_th;

    printf("[modelComputeCSlm] Computing coefficients in parallel\n");

    Clock *clock = Clock_new();
    Clock_tic(clock);

    // A single iteration of this loop can be abstracted away as:
    // computeCSPair(l, m, P_lm_th, precomp)

#pragma omp parallel
{
    if (omp_get_thread_num() == 1) printf("[modelComputeCSlmPrecompOMP] computing CSlm values using %d threads\n", omp_get_num_threads());
}


    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            computeCSPairOMP(l, m, data, precomp, C_lm, S_lm, P_lm_th);
        }
    }

    Clock_toc(clock);
    double time = elapsed_time(clock);
    free(clock);

    return time;

}

double modelComputeCSlmPrecompOMP2(SphericalModel *model, const data_iso *data, const Precomp *precomp) {

    // const int lmax = model->lmax;
    const int ll = model->ll;

    Matrix_d *C_lm = model->C_lm, *S_lm = model->S_lm;
    // Matrix_d *P_lm_th = precomp->Plm_th;

    printf("[modelComputeCSlm] Computing coefficients in parallel\n");

    Clock *clock = Clock_new();
    Clock_tic(clock);

    #pragma omp parallel
    {
        if (omp_get_thread_num() == 0) printf("[modelComputeCSlmPrecompOMP2] computing CSlm values using %d threads\n", omp_get_num_threads());

        #pragma omp for
        for (int i = 0; i < ll; i++) {
            lm_pair lm = i_to_lm(i);
            computeCSPair(lm.l, lm.m, data, precomp, C_lm, S_lm);
        }
    }

    Clock_toc(clock);
    double time = elapsed_time(clock);
    free(clock);

    return time;

}

double modelComputeCSlmPrecompOMP2Threads(SphericalModel *model, const data_iso *data, const Precomp *precomp, int nthreads) {

    // const int lmax = model->lmax;
    const int ll = model->ll;

    Matrix_d *C_lm = model->C_lm, *S_lm = model->S_lm;
    // Matrix_d *P_lm_th = precomp->Plm_th;

    printf("[modelComputeCSlm] Computing coefficients in parallel\n");
    omp_set_num_threads(nthreads);

    Clock *clock = Clock_new();
    Clock_tic(clock);

    #pragma omp parallel
    {
        if (omp_get_thread_num() == 0) printf("[modelComputeCSlmPrecompOMP2] computing CSlm values using %d threads\n", omp_get_num_threads());

        #pragma omp for
        for (int i = 0; i < ll; i++) {
            lm_pair lm = i_to_lm(i);
            computeCSPair(lm.l, lm.m, data, precomp, C_lm, S_lm);
        }
    }

    Clock_toc(clock);
    double time = elapsed_time(clock);
    free(clock);

    return time;

}

// Compute the value of Clm and Slm that are assigned to this process.
double modelComputeCSlmPrecompMPI(SphericalModel *model, const data_iso *data, const Precomp *precomp, int world_size, int this_rank) {

    // Total work is LL
    // const int lmax = model->lmax;
    const int ll = model->ll;

    // Only the ROOT RANK will actually fill the spherical model
    // Every other rank will, however, compute their own stretch of 

    Clock *clock = Clock_new();
    Clock_tic(clock);

    // Go ahead and compute the workload needed by this rank
    int this_workload = compute_workload(ll, world_size, this_rank);


    Matrix_i *start_end = compute_startend_array(ll, world_size);

    int i_start = vecat_i(start_end, this_rank) + 1;
    int i_end = vecat_i(start_end, this_rank + 1);

    // Allocate the space for Clm and Slm matrices
    Matrix_d *Clm = Matrix_new_d(1, this_workload);
    Matrix_d *Slm = Matrix_new_d(1, this_workload);

    // for (int i = 0; i < world_size; i++) {
    //     if (i == this_rank) {
    //         printf("[%d] workload: %d\n", this_rank, this_workload);
    //         printf("[%d] Allocated matrix of size %d\n", this_rank, this_workload);
    //         printf("[%d] Processing i \\in [%d, %d]\n", this_rank, i_start, i_end);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // Here we need to assert the validity of data and precomp
    assert(data != NULL);
    assert(precomp != NULL);
    assert(Clm != NULL);
    assert(Slm != NULL);

    // for (int i = 0; i < world_size; i++) {
    //     if (this_rank == i) {
    //         head_data(data);
    //     }

    //     MPI_Barrier(MPI_COMM_WORLD);
    // }   

    // Data looks legit, how about the precomp data? I'm worried that somethings getting fucked up there.
    // precomp data is legit...

    // printf("[%d] verified  data and precomp\n", this_rank);

    MPI_Barrier(MPI_COMM_WORLD);

    // Something isn't properly initialized before entering computeCSPair
    for (int i = i_start; i <= i_end; i++) {
        lm_pair pair = i_to_lm(i);
        // printf("computing CS(%d, %d)\n", pair.l, pair.m);
        computeCSPairMPI(pair.l, pair.m, data, precomp, Clm, Slm, i_start);
    }


    Matrix_i *recvcounts = compute_workload_array(ll, world_size);
    Matrix_i *displacements = compute_displacements(ll, world_size);

    // MPI_Barrier(MPI_COMM_WORLD);

    // for (int i = 0; i < world_size; i++) {
    //     if (this_rank == i) {
    //         printf("Clm: i \\in [%d, %d]\t", i_start, i_end);
    //         Matrix_print_d(Clm);
    //     }

    //     MPI_Barrier(MPI_COMM_WORLD);
    // }   

    // MPI_Barrier(MPI_COMM_WORLD); 

    // for (int i = 0; i < world_size; i++) {
    //     if (this_rank == i) {
    //         printf("Slm: i \\in [%d, %d]\t", i_start, i_end);
    //         Matrix_print_d(Slm);
    //     }

    //     MPI_Barrier(MPI_COMM_WORLD);
    // }    

    // MPI_Barrier(MPI_COMM_WORLD); 

    // Matrix_free_d(Clm);
    // Matrix_free_d(Slm);

    // MPI_Finalize();
    // exit(0);
    // Now we need to gatherv the results stored in Clm and Slm and add them spherical model

    // C coefficients
    MPI_Gatherv(Clm->data, 
        this_workload, 
        MPI_DOUBLE, 
        model->C_lm->data, 
        recvcounts->data, 
        displacements->data, 
        MPI_DOUBLE, 
        0, 
        MPI_COMM_WORLD);

    // S coefficients
    MPI_Gatherv(Slm->data, 
        this_workload, 
        MPI_DOUBLE, 
        model->S_lm->data, 
        recvcounts->data, 
        displacements->data, 
        MPI_DOUBLE, 
        0, 
        MPI_COMM_WORLD);

    Clock_toc(clock);
    
    double runtime = elapsed_time(clock);

    Matrix_free_i(recvcounts);
    Matrix_free_i(displacements);
    Matrix_free_d(Clm);
    Matrix_free_d(Slm);
    Matrix_free_i(start_end);
    free(clock);

    return runtime;
    // Hopefully this shit works lmaoo
}




// Write the coefficients of a spherical model to CSV
// Use the prefix to indicate what data set the spherical model was trained on
// pass "" as type if you don't want to specify
void SphericalModelToTXT(const SphericalModel *model, const char *type) {

    // Open a csv file
    char txt_file[100] = {0};

    if (strlen(type) > 0)
        sprintf(txt_file, "sph_%s_%d.txt", type, (int) model->lmax);
    else 
        sprintf(txt_file, "sph_%d.txt", (int) model->lmax);


    FILE *txt_out = fopen(txt_file, "w");

    for (size_t l = 0; l <= model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {
            fprintf(txt_out, "%lu\t%lu\t%lf\t%lf\n", l, m, model->C_lm->data[PT(l, m)], model->S_lm->data[PT(l, m)]);
        }
    }

    printf("[SphericalModelToTXT] wrote model to %s\n", txt_file);

    fclose(txt_out);

}

// Write the coefficients of a spherical model to CSV
// Use the prefix to indicate what data set the spherical model was trained on
// pass "" as type if you don't want to specify
// Binary format: LL * sizeof(double) bytes dedicated to the first LL C coefficients,
// then           LL * sizeof(double) bytes dedicated to the first LL S coefficients 
void SphericalModelToBIN(const SphericalModel *model, const char *type) {

    // Open a csv file
    char bin_file[100] = {0};

    if (strlen(type) > 0)
        sprintf(bin_file, "sph_%s_%d.bin", type, (int) model->lmax);
    else 
        sprintf(bin_file, "sph_%d.bin", (int) model->lmax);

    FILE *bin_out = fopen(bin_file, "wb");

    const size_t n_bytes = sizeof(double) * model->ll;

    fwrite(model->C_lm->data, n_bytes, 1, bin_out);
    fwrite(model->S_lm->data, n_bytes, 1, bin_out);

    fclose(bin_out);

}

SphericalModel *loadSphericalModel(const char *bin_in, int lmax) {

    if (lmax < 0) errx(1, "lmax [%d] must be positive\n", lmax);

    // Check to make sure file exists
    FILE *bin = fopen(bin_in, "rb");
    if (bin == NULL) errx(1, "Binary file %s does not exist\n", bin_in);

    SphericalModel *model = (SphericalModel *) malloc(sizeof(*model));

    if (model == NULL) errx(1, "Error allocating SphericalModel\n");

    // Allocate space for model
    model->lmax = lmax;
    model->ll = (lmax + 1) * (lmax + 2) / 2;
    model->C_lm = Matrix_new_d(1, model->ll);
    model->S_lm = Matrix_new_d(1, model->ll);

    const size_t n_bytes = sizeof(double) * model->ll;

    fread(model->C_lm->data, n_bytes, 1, bin);
    fread(model->S_lm->data, n_bytes, 1, bin);

    fclose(bin);

    return model;
}


void freeSphericalModel(SphericalModel *model) {

    if (model == NULL) return;

    if (model->C_lm) Matrix_free_d(model->C_lm);
    if (model->S_lm) Matrix_free_d(model->S_lm);

    free(model);

}