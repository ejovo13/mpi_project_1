/**========================================================================
 * ?                          model.c
 * @brief   : Functions dealing with the computational aspect of 
 *            _generating_ new models 
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-15
 *========================================================================**/

#include "geodesy.h"

spherical_model *copy_spherical_model(const spherical_model *rhs) {
    spherical_model *out = (spherical_model *) malloc(sizeof(*out));
    out->lmax    = rhs->lmax;
    out->ll      = rhs->ll;
    out->C_lm    = matclone_d(rhs->C_lm);
    out->S_lm    = matclone_d(rhs->S_lm);
    out->P_lm_th = matclone_d(rhs->P_lm_th);
    out->pcs     = matclone_d(rhs->pcs);
    out->clm     = matclone_d(rhs->clm);
    out->slm     = matclone_d(rhs->slm);
    return out;
}

iso_model *copy_iso (const iso_model *rhs) {
    iso_model *out = (iso_model *) malloc(sizeof(*out));
    out->data = rhs->data; // underlying data DOES NOT CHANGE
    out->coeff = rhs->coeff; // underlying coefficients doesn't change either
    out->model = copy_spherical_model(rhs->model);
    return out;
};

// copy the current state of this iso, deep copying all of its components
iso_trained *copy_trained_model (const iso_trained *rhs) {

    iso_trained *out = (iso_trained *) malloc(sizeof(*out));
    out->Clm_history = matclone_d(rhs->Clm_history);
    out->Slm_history = matclone_d(rhs->Slm_history);
    out->MSE         = matclone_d(rhs->MSE);
    out->dMSE        = matclone_d(rhs->dMSE);
    out->count       = rhs->count;
    out->iso         = copy_iso(out->iso);
    return out;
};

void free_spherical_model(spherical_model *sph_model) {

    if (sph_model == NULL) return;

    if (sph_model->C_lm != NULL)
        Matrix_free_d(sph_model->C_lm);

    if (sph_model->S_lm != NULL)
        Matrix_free_d(sph_model->S_lm);

    if (sph_model->P_lm_th != NULL)
        Matrix_free_d(sph_model->P_lm_th);

    if (sph_model->pcs != NULL)
        Matrix_free_d(sph_model->pcs);

    if (sph_model->clm != NULL)
        Matrix_free_d(sph_model->clm);

    if (sph_model->slm != NULL)
        Matrix_free_d(sph_model->slm);

}

iso_model *newModel(data_iso *data, int lmax) {

    if (lmax < 0) {
        perror("Negative lmax passed into newModel, exiting\n");
        exit(EXIT_FAILURE);
    }

    const int ll = (lmax + 1) * (lmax + 2) / 2;

    spherical_model *model = (spherical_model *) malloc(sizeof(*model));

    if (model == NULL) {
        perror("Error allocating new model\n");
        exit(EXIT_FAILURE);
    }

    model->C_lm    = NULL;
    model->S_lm    = NULL;
    model->P_lm_th = NULL;
    model->pcs     = NULL;
    model->clm     = NULL;
    model->slm     = NULL;

    model->ll = ll;
    model->lmax = lmax;
    model->C_lm    = Matrix_new_d(1, model->ll);
    model->S_lm    = Matrix_new_d(1, model->ll);
    model->P_lm_th = Matrix_new_d(data->t, model->ll);

    // Now process the data to fill in the model    
    // Initialize P_lm_th

    // printf("th: [%lf, %lf]\tph: [%lf, %lf]\n", th_0, th_f, ph_0, ph_f);

    // make sure data got loaded
    printf("[newModel] New model successfully constructed\n");

    printf("[newModel] Model: { .ll = %lu, .lmax = %lu}\n", model->ll, model->lmax);

    printf("[newModel] Data:  { .N = %d, .t = %d, .p = %d}\n", data->N, data->t, data->p);

    // Now let's create the "iso_model type"
    iso_model *iso = (iso_model *) malloc(sizeof(*iso));

    iso->data = data;
    iso->model = model;

    modelComputePlm(iso, data);
    printf("[newModel] Completed computing Plm\n");
    modelComputeCSlm(model, data);
    printf("[newModel] Completed computing CSlm\n");

    // initialize coefficients
    iso->model->clm = compute_mse_coeff_clm(data, model); // 2 x ll matrix
    iso->model->slm  = compute_mse_coeff_slm(data, model); // 2 x ll matrix
    iso->model->pcs = compute_pcs(data, model, iso->model->clm, iso->model->slm); // sum(clm * Plmcos * cos)_i for i in 1..N

    return iso;

}

// Compute a new model, storing the time that it took at the address of `time_taken`
iso_model *newModelTimed(data_iso *data, int lmax, double *time_taken) {

    if (lmax < 0) {
        perror("Negative lmax passed into newModel, exiting\n");
        exit(EXIT_FAILURE);
    }

    const int ll = (lmax + 1) * (lmax + 2) / 2;

    spherical_model *model = (spherical_model *) malloc(sizeof(*model));

    if (model == NULL) {
        perror("Error allocating new model\n");
        exit(EXIT_FAILURE);
    }

    model->C_lm    = NULL;
    model->S_lm    = NULL;
    model->P_lm_th = NULL;
    model->pcs     = NULL;
    model->clm     = NULL;
    model->slm     = NULL;

    model->ll = ll;
    model->lmax = lmax;
    model->C_lm    = Matrix_new_d(1, model->ll);
    model->S_lm    = Matrix_new_d(1, model->ll);
    model->P_lm_th = Matrix_new_d(data->t, model->ll);

    // Now process the data to fill in the model    
    // Initialize P_lm_th

    // printf("th: [%lf, %lf]\tph: [%lf, %lf]\n", th_0, th_f, ph_0, ph_f);

    // make sure data got loaded
    printf("[newModel] New model successfully constructed\n");

    printf("[newModel] Model: { .ll = %lu, .lmax = %lu}\n", model->ll, model->lmax);

    printf("[newModel] Data:  { .N = %d, .t = %d, .p = %d}\n", data->N, data->t, data->p);

    // Now let's create the "iso_model type"
    iso_model *iso = (iso_model *) malloc(sizeof(*iso));

    iso->data = data;
    iso->model = model;

    modelComputePlm(iso, data);
    printf("[newModel] Completed computing Plm\n");
    *time_taken = modelComputeCSlm(model, data);
    printf("[newModel] Completed computing CSlm\n");

    // initialize coefficients
    iso->model->clm = compute_mse_coeff_clm(data, model); // 2 x ll matrix
    iso->model->slm  = compute_mse_coeff_slm(data, model); // 2 x ll matrix
    iso->model->pcs = compute_pcs(data, model, iso->model->clm, iso->model->slm); // sum(clm * Plmcos * cos)_i for i in 1..N
    iso->f_hat = NULL;

    return iso;

}


void writeModel(const spherical_model *model, const data_iso *data, const char *prefix) {

    // const int ll = model->ll;

    char fileout[100] = {0};

    sprintf(fileout, "%smodel_%lu_%d_%d.txt", prefix, model->lmax, data->t, data->p);

    //model_lmax_nth_nph.txt
    printf("\nWriting data to file: %s\n\n", fileout);

    FILE *out = fopen(fileout, "w");

    for (size_t l = 0; l <= model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {
            fprintf(out, "%lu\t%lu\t%.15lf\t%.15lf\n", l, m, model->C_lm->data[PT(l, m)], model->S_lm->data[PT(l, m)]);
        }
    }

    fclose(out);

}

// Write an artifically modified model where the written Clm and Slm are multiplied by 4 if 
// m != 0
void writeModifiedModel(const spherical_model *model, const data_iso *data, const char *prefix) {

    // const int ll = model->ll;

    char fileout[100] = {0};

    sprintf(fileout, "%smmodel_%lu_%d_%d.txt", prefix, model->lmax, data->t, data->p);

    //model_lmax_nth_nph.txt
    printf("\nWriting data to file: %s\n\n", fileout);

    FILE *out = fopen(fileout, "w");

    for (size_t l = 0; l <= model->lmax; l++) {
        fprintf(out, "%lu\t%u\t%.15lf\t%.15lf\n", l, 0, model->C_lm->data[PT(l, 0)], model->S_lm->data[PT(l, 0)]);
        for (size_t m = 1; m <= l; m++) {
            fprintf(out, "%lu\t%lu\t%.15lf\t%.15lf\n", l, m, 4.0 * model->C_lm->data[PT(l, m)], 4.0 * model->S_lm->data[PT(l, m)]);
        }
    }

    fclose(out);

}

void modelComputePlm(iso_model *iso, const data_iso *data) {

    spherical_model *model = iso->model;

    // const int ll = model->ll;
    Matrix_d *P_lm_th = Matrix_new_d(model->P_lm_th->nrows, model->P_lm_th->ncols);
    Matrix_free_d(model->P_lm_th);
    model->P_lm_th = P_lm_th;
    

    printf("[modelComputePlm] Setting up spherical harmonics with lmax: %lu\n", model->lmax);

    struct spherical_harmonics *sph_model = (struct spherical_harmonics *) malloc(sizeof(*sph_model)); // used strictly to compute the P_lm_th matrix
	setup_spherical_harmonics(model->lmax, sph_model);

    for (int i = 0; i < data->t; i++) { // store this data in a more compact matrix

        // P is a big matrix who stores the information as (p00(th0) p10(th0) )
        // printf("Processing i: %d, ptr: %p\n", i, matacc_d(model->P_lm_th, i, 0));
		computeP(sph_model, matacc_d(model->P_lm_th, i, 0), cos(vecat_d(data->th, i))); // should be mathematically equal to P(sin(phi_p)), angle named

        // // double sinth = sin(data->th[i]);
        // double sinth = sin(vecat_d(data->th, i));
        // for (int j = 0; j < ll; j++) {
        //     *matacc_d(model->P_lm_th, i, j) *= sinth;
        // }
    }

    iso->coeff = sph_model;    
    // Matrix_print_d(model->P_lm_th);
}

// Estimate the Laplace series coefficients Clm and Slm via numeric integration,
// and then return the time that it took to compute this value.
double modelComputeCSlm(spherical_model *model, const data_iso *data) {

    const int lmax = model->lmax;
    // const int ll = model->ll;

    Matrix_d *C_lm = model->C_lm, *S_lm = model->S_lm;
    // Matrix_d *P_lm = model->P_lm_th;

    printf("[ modelComputeCSlm ] :\n");
    // Matrix_print_d(P_lm);


    // Initialize integral values
    double c_integral = 0;
    double s_integral = 0;
    int count = 0;

    Clock *clock = Clock_new();
    Clock_tic(clock);

    for (int l = 0; l <= lmax; l++) {

        for (int m = 0; m <= l; m++) {

            // Now that I have the Associated legendre functions and the data efficiently loaded, let's write
            // the code to approximate the integrals
            c_integral = 0;
            s_integral = 0;

            // compute vector of sinth to reduce computational workload
            Matrix_d *sinth = Matrix_new_d(1, data->t);
            for (int i_th = 0; i_th < data->t; i_th++) {
                *vecptr_d(sinth, i_th) = sin(vecat_d(data->th, i_th));
            }

            for (int i = 0; i < data->N; i++) {

                int i_th = i % data->t;
                int i_ph = i / data->t;

                // if (l == 0 && m == 0) {
                    // printf("(i = %d) Calculating (th, ph) : f(%lf, %lf) = %lf\t", i, vecat_d(data->th, i_th), vecat_d(data->ph, i_ph), matat_d(data->r, i_ph, i_th)) ;
                    // printf("= f(%lf, %lf) [f(phi, lambda)]\n", thew_to_phip(vecat_d(data->th, i_th)), 
                                                            //    phiw_to_lamp(vecat_d(data->ph, i_ph)));
                // }

                double ph_iph = vecat_d(data->ph, i_ph);
                double cos_mph = cos(m * ph_iph); // could use a recursive relationship to calclute this faster
                double sin_mph = sin(m * ph_iph);

                // use the midpoint formula
                c_integral += data->r[i] * matat_d(model->P_lm_th, i_th, PT(l, m)) * cos_mph * vecat_d(sinth, i_th);
                s_integral += data->r[i] * matat_d(model->P_lm_th, i_th, PT(l, m)) * sin_mph * vecat_d(sinth, i_th);

                count ++;
            }

            c_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);
            s_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);

            // c_integral *= 1.023491;
            // s_integral *= 1.023491;
            // c_integral *= (data->dp * data->dt) / sqrt((2.0 * TWO_PI));
            // s_integral *= (data->dp * data->dt) / sqrt((2.0 * TWO_PI));
            // c_integral *= (data->dp * data->dt) / (4.0 * TWO_PI);
            // s_integral *= (data->dp * data->dt) / (4.0 * TWO_PI);
            // c_integral *= data->dp * data->dt / (sqrt(TWO_PI));
            // s_integral *= data->dp * data->dt / (sqrt(TWO_PI));

            // But I actually think that this shit is already normalized
            vecset_d(C_lm, PT(l, m), c_integral);
            vecset_d(S_lm, PT(l, m), s_integral);

            // printf("Computed coefficients for (%d, %d)\n", l, m);
        }
    }

    Clock_toc(clock);
    double time = elapsed_time(clock);
    free(clock);

    return time;

}


// Given an lmax, compute the values of C_lm and S_lm and return the values
// in a Vector of the form
//
//  Coeff = [ C_lm S_lm ]
//
// where the length of Coeff is 2 * (lmax + 1) * (lmax + 2) / 2
iso_model *compute_model(int lmax, const char *datafile, int npoints) {

    const int ll = (lmax + 1) * (lmax + 2) / 2;
    const int n_theta = sqrt(npoints / 2);
    const int n_phi   = 2 * n_theta;

    printf("[compute_model] Constructing model { lmax: %d, t: %d, p: %d}\n", lmax, n_theta, n_phi);
    data_iso *data = load_data_iso(datafile, n_theta, n_phi);

    // Now let's print the data points that were actually loaded.

    // TODO refactor to be stored in the model
    // use lmax to initialize data->l_indices 
    Matrix_i *l_indices = Matrix_new_i(1, ll);
    Matrix_i *m_indices = Matrix_new_i(1, ll);

    int i = 0;
    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            l_indices->data[i] = l;
            m_indices->data[i] = m;
            i++;
        }
    }

    iso_model *iso = newModel(data, lmax);
    writeModel(iso->model, data, "");    

    return iso;
}

// Given an lmax, compute the values of C_lm and S_lm and return the values
// in a Vector of the form
//
//  Coeff = [ C_lm S_lm ]
//
// where the length of Coeff is 2 * (lmax + 1) * (lmax + 2) / 2
iso_model *compute_model_binary(int lmax, const char *binary_in, int npoints, bool __log) {

    const int ll = (lmax + 1) * (lmax + 2) / 2;
    const int n_theta = sqrt(npoints / 2);
    const int n_phi   = 2 * n_theta;

    if (__log) 
        printf("[compute_model] Constructing model { lmax: %d, t: %d, p: %d}\n", lmax, n_theta, n_phi);
    data_iso *data = load_data_binary(binary_in, n_theta, n_phi, __log);

    // Now let's print the data points that were actually loaded.


    // use lmax to initialize data->l_indices 
    Matrix_i *l_indices = Matrix_new_i(1, ll);
    Matrix_i *m_indices = Matrix_new_i(1, ll);

    int i = 0;
    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            l_indices->data[i] = l;
            m_indices->data[i] = m;
            i++;
        }
    }

    iso_model *iso = newModel(data, lmax);
    writeModel(iso->model, data, "");    

    return iso;
}

// I'm not particularly fond of this inconsistent naming scheme between
// pascalCase and snake_case
void modelFree(iso_model* iso) {

    if (iso == NULL) return;

    free_data_iso(iso->data);
    free_spherical_model(iso->model);
    free_spherical_harmonics(iso->coeff);

    if (iso->f_hat != NULL) 
        Matrix_free_d(iso->f_hat);

    free(iso);
}

// Return the time that it takes to compute a model of size Lmax using 
// numerical quadrature
double time_new_model(int lmax, int npoints, const char *data_filename) {

    double time = 0;

    // First step is to load the model
    const int ll = (lmax + 1) * (lmax + 2) / 2;
    const int n_theta = sqrt(npoints / 2);
    const int n_phi   = 2 * n_theta;

    printf("[compute_model] Constructing model { lmax: %d, t: %d, p: %d}\n", lmax, n_theta, n_phi);
    data_iso *data = load_data_iso(data_filename, n_theta, n_phi);

    // Now let's print the data points that were actually loaded.


    // use lmax to initialize data->l_indices 
    Matrix_i *l_indices = Matrix_new_i(1, ll);
    Matrix_i *m_indices = Matrix_new_i(1, ll);

    int i = 0;
    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            l_indices->data[i] = l;
            m_indices->data[i] = m;
            i++;
        }
    }

    iso_model *iso = newModelTimed(data, lmax, &time);
    writeModel(iso->model, data, "timed");    

    // free the model...
    modelFree(iso); 

    return time;
}



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
Matrix_d *compute_mse_coeff(const data_iso *data, const spherical_model *model) {

    // const int lmax = Matrix_size_d(P_lm_th);
    // const int 
    // const int ll = (lmax + 1) * (lmax + 2) / 2;
    const int ll = model->ll;
    const int lmax = model->lmax;



    // double *CS = (double *) malloc(sizeof(*CS) * 2 * ll * data->N); // start off using full sized array
    Matrix_d *cs_mat = Matrix_new_d(data->N * 2, ll);
    // double *C = cs_mat->data;
    // double *S = matacc_d(cs_mat, 1, 0); // halfway point

    // int index = 0;;
    double cosmph = 0;
    double sinmph = 0;
    const int N = data->N;

    // for (int i = 0; i < data->N; i++) {
    for (int i = 0; i < data->N; i++) {

        for (int l = 0; l <= lmax; l++) {
            // th changes every single cycle
            for (int m = 0; m <= l; m++) {

                int j = PT(l, m);
                 
                int i_ph = i / data->t; // ph has floor division cycles
                int i_th = i % data->t;
                
                cosmph = cos(m * vecat_d(data->ph, i_ph));
                sinmph = sin(m * vecat_d(data->ph, i_ph)); 

                *matacc_d(cs_mat, i, j)     = matat_d(model->P_lm_th, i_ph, i_th) * cosmph;
                *matacc_d(cs_mat, i + N, j) = matat_d(model->P_lm_th, i_ph, i_th) * sinmph;

                // printf("(%d) Computed mse coefficients c(%d, %d)\n", i, l, m);
            }
        }
        // printf("(i = %d) clm_i = ", i);
        // MatIter_print_d(Matrix_row_begin_d(cs_mat, i), Matrix_row_end_d(cs_mat, i));
        
    }

    return cs_mat;
}

Matrix_d *compute_mse_coeff_clm(const data_iso *data, const spherical_model* model) {

    Matrix_d *clm = Matrix_new_d(data->N, model->ll);
    double cosmph = 0;

    // printf("lmax: %lu\n", model->lmax);

    for (int i = 0; i < data->N; i++) {
        for (size_t l = 0; l <= model->lmax; l++) {
            // th changes every single cycle
            for (size_t m = 0; m <= l; m++) {

                int j = PT(l, m);
                 
                int i_ph = i / data->t; // ph has floor division cycles
                int i_th = i % data->t;

                // printf("j: %d\n", j);
                
                // cosmph = cos(m * data->ph[i_ph]); 
                cosmph = cos(m * vecat_d(data->ph, i_ph)); 
                *matacc_d(clm, i, j) = matat_d(model->P_lm_th, i_th, j) * cosmph;
            }
        }
    }

    return clm;
}

Matrix_d *compute_mse_coeff_slm(const data_iso *data, const spherical_model* model) {

    Matrix_d *slm = Matrix_new_d(data->N, model->ll);
    double sinmph = 0;

    for (int i = 0; i < data->N; i++) {
        for (size_t l = 0; l <= model->lmax; l++) {
            // th changes every single cycle
            for (size_t m = 0; m <= l; m++) {

                int j = PT(l, m);
                 
                int i_ph = i / data->t; // ph has floor division cycles
                int i_th = i % data->t;
                
                // sinmph = cos(m * data->ph[i_ph]); 
                sinmph = sin(m * vecat_d(data->ph, i_ph)); 
                *matacc_d(slm, i, j)     = matat_d(model->P_lm_th, i_th, j) * sinmph;
            }
        }
    }

    return slm;
}


// O(lmax^2)
double compute_pcs_i(const spherical_model *model, const Matrix_d *clm, const Matrix_d *slm, int i) {

    // I already have the values of c_li^m and s_li^m, just need to
    // combine them with the parameters of my model
    double sum = 0;

    for (size_t l = 0; l <= model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {

            int index = PT(l, m);

            sum += model_C(model, l, m) * matat_d(clm, i, index);
            sum += model_S(model, l, m) * matat_d(slm, i, index);
        }
    }

    return sum;
}

Matrix_d *compute_pcs(const data_iso *data, const spherical_model *model, const Matrix_d *clm, const Matrix_d *slm) {

    // Pcs will simply have length data->N
    // double * pcs = (double *) malloc(sizeof(*pcs) * data->N);
    Matrix_d *pcs_mat = Matrix_new_d(1, data->N);

    // for i in {1..N} compute PCS_i 
    for (int i = 0; i < data->N; i++) {
        pcs_mat->data[i] = compute_pcs_i(model, clm, slm, i);
    }

    return pcs_mat;
}

Matrix_d *compute_pcs_indices(const iso_model *iso, const Matrix_i *indices) {

    const int n = Matrix_size_i(indices);
    // Pcs will simply have length data->N
    Matrix_d *pcs_mat = Matrix_new_d(1, n);

    // for i in {1..N} compute PCS_i 
    for (int __i = 0; __i < n; __i++) {
        int i = indices->data[__i];
        pcs_mat->data[__i] = compute_pcs_i(iso->model, iso->model->clm, iso->model->slm, i);
    }

    return pcs_mat;
}

// Compute the gradient
double compute_gradient_clm(const data_iso *data, const Matrix_d *clm, const Matrix_d *pcs, int l, int m) {

    double sum = 0;

    for (int i = 0; i < data->N; i++) {

        double clmi = matat_d(clm, i, PT(l, m));
        sum += 2 * clmi * pcs->data[i] - 2 * clmi * data->r[i];
    }

    return sum;
}

// __i corresponds to a simple counter whereas i corresponds to the __ith index in `indices`
double compute_gradient_clm_points(const iso_model *iso, const Matrix_d *clm, const Matrix_d *pcs, int l, int m, const Matrix_i *indices) {

    double sum = 0;

    // iterate through the indices
    for (size_t __i = 0; __i < Matrix_size_i(indices); __i++) {
        int i = indices->data[__i];
        double clmi = matat_d(clm, i, PT(l, m));
        sum += 2 * clmi * pcs->data[__i] - 2 * clmi * iso->data->r[i];
    }

    return sum;
}

// Passed pcs is a REDUCED estimation of PCS
double compute_gradient_slm_points(const iso_model *iso, const Matrix_d *slm, const Matrix_d *pcs, int l, int m, const Matrix_i *indices) {

    double sum = 0;

    // iterate through the indices
    for (size_t i = 0; i < Matrix_size_i(indices); i++) {
        int index = indices->data[i];
        double slmi = matat_d(slm, index, PT(l, m));
        sum += 2 * slmi * pcs->data[i] - 2 * slmi * iso->data->r[index];
    }

    return sum;
}

double compute_gradient_slm(const data_iso *data, const Matrix_d *slm, const Matrix_d *pcs, int l, int m) {


    double sum = 0;

    for (int i = 0; i < data->N; i++) {

        double slmi = matat_d(slm, i, PT(l, m));
        // sum += 2 * slmi * pcs->data[i] - 2 * slmi * data->r[i];
        sum += 2 * slmi * pcs->data[i] - 2 * slmi * data->r[i];
    }

    return sum;
}

