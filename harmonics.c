#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <assert.h>
#include <sys/time.h>
#include <inttypes.h>

#include "harmonics.h"

double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, NULL);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

/* represent n in <= 8 char  */
void human_format(char * target, long n) {
        if (n < 1000) {
                sprintf(target, "%" PRId64, n);
                return;
        }
        if (n < 1000000) {
                sprintf(target, "%.1f K", n / 1e3);
                return;
        }
        if (n < 1000000000) {
                sprintf(target, "%.1f M", n / 1e6);
                return;
        }
        if (n < 1000000000000ll) {
                sprintf(target, "%.1f G", n / 1e9);
                return;
        }
        if (n < 1000000000000000ll) {
                sprintf(target, "%.1f T", n / 1e12);
                return;
        }
}


void setup_spherical_harmonics(int lmax, struct spherical_harmonics *self)
{
	self->lmax = lmax;
	int sizeCS = (lmax + 1) * (lmax + 1);
	int sizeAB = (lmax + 1) * (lmax + 2) / 2;

	self->CS = malloc(sizeCS * sizeof(double));
	self->A = malloc(sizeAB * sizeof(double));
	self->B = malloc(sizeAB * sizeof(double));
	if (self->CS == NULL || self->A == NULL || self->B == NULL)
		err(1, "cannot allocate harmonics\n");

	/* compute the A, B coefficients */
	for (int l = 2; l <= lmax; l++) {
		double ls = l * l;
		double lm1s = (l - 1) * (l - 1);
		for (int m = 0; m < l - 1; m++) {
			double ms = m * m;
			self->A[PT(l, m)] = sqrt((4 * ls - 1) / (ls - ms));
			self->B[PT(l, m)] = -sqrt((lm1s - ms) / (4 * lm1s - 1));
            // printf("P(%d, %d) -> %lf ", l, m, )
		}
	}

    // print the computed values:


}

void load_data_points(const char *filename, int npoint, struct data_points *self)
{
	self->npoint = npoint;
	self->phi = malloc(npoint * sizeof(double));
	self->lambda = malloc(npoint * sizeof(double));
	self->V = malloc(npoint * sizeof(double));
	if (self->phi == NULL || self->lambda == NULL || self->V == NULL)
		err(1, "cannot allocate data points\n");

	FILE *f = fopen(filename, "r");
	if (f == NULL)
		err(1, "cannot open %s", filename);
	
	for (int i = 0; i < npoint; i++) {
		int k = fscanf(f, "%lg %lg %lg", &self->lambda[i], &self->phi[i], &self->V[i]);
		if (k == EOF) {
			if (ferror(f))
				err(1, "read error");
			errx(1, "premature end-of-file after %d aecords", i);
		}
		if (k != 3)
			errx(1, "parse error on line %d", i+1);
	}
	fclose(f);
}

data_points_iso *load_data_points_iso(const char *filename, int npoint)
{

    data_points_iso *data = (data_points_iso *) malloc(sizeof(*data));

    if (!data) err(1, "Cannot allocate data points structure");

	data->npoint = npoint;
	data->th = malloc(npoint * sizeof(double));
	data->ph = malloc(npoint * sizeof(double));
	data->r  = malloc(npoint * sizeof(double));
	if (data->th == NULL || data->ph == NULL || data->r == NULL)
		err(1, "cannot allocate data points\n");

	FILE *f = fopen(filename, "r");
	if (f == NULL)
		err(1, "cannot open %s", filename);

    double p_lambda;
    double p_phi;
    double p_value;    

	for (int i = 0; i < npoint; i++) {

		// int k = fscanf(f, "%lg %lg %lg", &data->lambda[i], &data->phi[i], &data->V[i]);
		int k = fscanf(f, "%lg %lg %lg", &p_lambda, &p_phi, &p_value);

        // Transformation from NASA's format to ISO format
        if (p_lambda < 0) {
            data->ph[i] = p_lambda + TWO_PI;
        } else {
            data->ph[i] = p_lambda;
        }
        data->th[i] = HALF_PI - p_phi;
        data->r[i] = p_value;

		if (k == EOF) {
			if (ferror(f))
				err(1, "read error");
			errx(1, "premature end-of-file after %d records", i);
		}
		if (k != 3)
			errx(1, "parse error on line %d", i+1);
	}
	fclose(f);

    return data;
}

void print_npoints(const data_points_iso* data, int npoints) {

    printf("th\t\tph\t\tr\n");

    for (int i = 0; i < npoints; i++) {
        printf("%lf\t%lf\t%lf\n", data->th[i], data->ph[i], data->r[i]);
    }

}

void write_npoints(const data_points_iso* data, int npoints, const char *filename) {

    FILE *out = fopen(filename, "w");

    for (int i = 0; i < npoints; i++) {
        fprintf(out, "%lf\t%lf\t%lf\n", data->th[i], data->ph[i], data->r[i]);
    }

    fclose(out);

}

// reproduce iso.csv using a limited amount of space
void write_iso(const data_iso* data, const char *filename) {

    FILE *out = fopen(filename, "w");

    for (int i = 0; i < data->N; i++) {
        // fprintf(out, "%lf\t%lf\t%lf\n", data->th[i % data->t], data->ph[i / data->p], data->r[i]);
        fprintf(out, "%.15lf\t%.15lf\t%.15lf\n", vecat_d(data->th, i % data->t), vecat_d(data->ph, i / data->t), vecat_d(data->r, i));
    }

    fclose(out);

}

void load_spherical_harmonics(const char *filename, int lmax, struct spherical_harmonics *self)
{
	FILE *g = fopen(filename, "r");

    int count = 0;

	if (g == NULL)
		err(1, "cannot open %s", filename);
	setup_spherical_harmonics(lmax, self);
	for (;;) {
		int l, m;
		double c, s;
		int k = fscanf(g, "%d %d %lg %lg", &l, &m, &c, &s);
		if (m == 0 && s != 0)
			errx(1, "non-zero S coefficient with l=%d and m=0", l);
		self->CS[CT(l, m)] = c;
		if (m > 0)
			self->CS[ST(l, m)] = s;
		if (k == EOF) {
			if (ferror(g))
				err(1, "read error");
			break;
		}
		if (k != 4) {
            printf("Error on line %d\n", count);
			errx(1, "parse error yabish");
        }

        // printf("Read: %d %d %lg %lg\n", l, m, c, s);

        count ++;
	}
	fclose(g);
}

// Load an iso model given Clm and Slm data points
// iso_model *

/*
 * Compute all the (fully normalized) Associated Legendre function of degree <= lmax.
 * On exit, P_{l,m} (with 0 <= m <= l <= lmax) can be found in P[PT(l, m)].
 * P must be preallocated of size (lmax * lmax + 3) / 2.
 * The constants A and B must be precomputed.
 */
void computeP(const struct spherical_harmonics *self, double *P, double sinphi)
{
	double cosphi = sqrt(1 - sinphi * sinphi);
	double temp = 1;
	P[PT(0, 0)] = temp;
	if (self->lmax == 0)
		return;
	P[PT(1, 0)] = sinphi * sqrt(3) * temp;
	temp = 0.5 * sqrt(3) * cosphi * temp;
	P[PT(1, 1)] = temp;
	for (int l = 2; l <= self->lmax; l++) {
		for (int m = 0; m < l - 1; m++)
			P[PT(l, m)] = self->A[PT(l, m)] * (sinphi * P[PT(l - 1, m)] + self->B[PT(l, m)] * P[PT(l - 2, m)]);
		P[PT(l, l - 1)] = sinphi * sqrt(2 * (l - 1) + 3) * temp;
		temp = -sqrt(1.0 + 0.5 / l) * cosphi * temp;
		P[PT(l, l)] = temp;
	}

    // for (int l = 0; l <= self->lmax; l++) {
    //     for (int m = 0; m <= l; m++) {
    //         printf("P(%d, %d) -> %lf ", l, m, P[PT(l, l - 1)]);
    //     }
    // }
}

// Predict any single arbitrary point f(\theta, \phi)
// Runs in O(l^2) time 
double modelPredict(const iso_model *iso, double theta, double phi) {

    // Compute the P values
    Matrix_d *P_lm = Matrix_new_d(1, iso->model->ll);

    // computeP(iso->coeff, P_lm->data, cos(theta));
    computeP(iso->coeff, P_lm->data, cos(theta));

    // Now that we have plm, compute the sum using the model's coefficients

    double sum = 0;

    for (size_t l = 0; l <= iso->model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {
            // Plm(cos\theta)*[clm cos(m\phi) + slm sin(m\phi)]
            sum += vecat_d(P_lm, PT(l, m)) * (model_C(iso->model, l, m) * cos(m * phi) +
                                              model_S(iso->model, l, m) * sin(m * phi));

            // printf("Summing Clm: %lf and Slm: %lf\n", model_C(iso->model, l, m), model_S(iso->model, l, m));
        }
    }

    // for (int i = 0; i < iso->model->ll; i++) {
        // sum += vecat_d(P_lm, i) * (vecat_d(iso->model->C_lm, i) * cos(m * phi))
    // }

    return sum;

}

// Predict a single point that belongs to the data set
double modelPredictDataPoint(const iso_model *iso, int i) {

    const Matrix_d *P_lm_th = iso->model->P_lm_th;

    int i_ph = data_i_ph(iso->data, i);
    int i_th = data_i_th(iso->data, i);

    double ph = vecat_d(iso->data->th, i_ph);
    double th = vecat_d(iso->data->ph, i_th);

    double sum = 0;

    for (size_t l = 0; l <= iso->model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {
            // Plm(cos\theta)*[clm cos(m\phi) + slm sin(m\phi)]
            sum += matat_d(P_lm_th, i_th, PT(l, m)) * (model_C(iso->model, l, m) * cos(m * ph)) +
                                              model_S(iso->model, l, m) * sin(m * ph);

            // printf("Summing Clm: %lf and Slm: %lf\n", model_C(iso->model, l, m), model_S(iso->model, l, m));
        }
    }

    return sum;

}

Matrix_d *modelPredictDataPoints(const iso_model *iso, Matrix_i *indices) {

    const int n = Matrix_size_i(indices);
    Matrix_d *predictions = Matrix_new_d(1, n);

    for (int i = 0; i < n; i++) {
        predictions->data[i] = modelPredictDataPoint(iso, indices->data[i]);
    }

    return predictions;
}

// Predict the values for a vector of indices passed in as argument
// Matrix_d *modelPredictDataPoints(const Matrix_i indices)

// Convert a data set to an output predicted by the model
Matrix_d *modelPredictData(const iso_model *iso) {

    Matrix_d *m = Matrix_new_d(iso->data->p, iso->data->t);

    // For all of the data points
    for (int i = 0; i < iso->data->N; i++) {
        m->data[i] = modelPredictDataPoint(iso, i);
    }

    return m;
}

// Compute the first n predictions of models data
Matrix_d *modelPredictN(const iso_model *iso, int __n) {

    int n = iso->data->N < __n ? iso->data->N : __n;

    const Matrix_d *th = iso->data->th;
    const Matrix_d *ph = iso->data->ph;

    Matrix_d *predictions = Matrix_new_d(1, n);

    for (int i = 0; i < n; i++) {

        int i_ph = i / iso->data->t;
        int i_th = i % iso->data->t;

        predictions->data[i] = modelPredict(iso, vecat_d(th, i_th), vecat_d(ph, i_ph));
        // printf("=======================\n");
    }

    return predictions;
}

double evaluate(const struct spherical_harmonics *self, const double *P, double lambda)
{
	int lmax = self->lmax;
    // int sizeCS = (lmax + 1) * (lmax + 2);
	int sizeCS = (lmax + 1) * (lmax + 1);
	double scratch[sizeCS];

    // Instead of using a weird scheme, just create two matrices one for C and one for S

	for (int l = 0; l <= lmax; l++) {
		/* zonal term */
		scratch[CT(l, 0)] = P[PT(l, 0)];

		/* tesseral terms */
		for (int m = 1; m <= l; m++) {
			scratch[CT(l, m)] = P[PT(l, m)] * cos(m * lambda);
			scratch[ST(l, m)] = P[PT(l, m)] * sin(m * lambda);
		}
	}

	/* dot product */
	double V = 0;	
	for (int i = 0; i < sizeCS; i++)
		V += scratch[i] * self->CS[i];
	return V;
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

// Create a uniform model so that I can compute some predictions
// iso_model *newModelUniform()

void writeModel(const spherical_model *model, const data_iso *data, const char *prefix) {

    const int ll = model->ll;

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

void modelComputePlm(iso_model *iso, const data_iso *data) {

    spherical_model *model = iso->model;

    const int ll = model->ll;
    Matrix_d *P_lm_th = Matrix_new_d(model->P_lm_th->nrows, model->P_lm_th->ncols);
    Matrix_free_d(model->P_lm_th);
    model->P_lm_th = P_lm_th;
    

    printf("[modelComputePlm] Setting up spherical harmonics with lmax: %lu\n", model->lmax);

    struct spherical_harmonics *sph_model = (struct spherical_harmonics *) malloc(sizeof(*sph_model)); // used strictly to compute the P_lm_th matrix
	setup_spherical_harmonics(model->lmax, sph_model);

    for (int i = 0; i < data->t; i++) { // store this data in a more compact matrix

        // P is a big matrix who stores the information as (p00(th0) p10(th0) )
        printf("Processing i: %d, ptr: %p\n", i, matacc_d(model->P_lm_th, i, 0));
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

// Estimate the Laplace series coefficients Clm and Slm via numeric integration
void modelComputeCSlm(spherical_model *model, const data_iso *data) {

    const int lmax = model->lmax;
    const int ll = model->ll;

    Matrix_d *C_lm = model->C_lm, *S_lm = model->S_lm;
    Matrix_d *P_lm = model->P_lm_th;

    printf("[ modelComputeCSlm ] :\n");
    // Matrix_print_d(P_lm);


    // Initialize integral values
    double c_integral = 0;
    double s_integral = 0;
    int count = 0;

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
                c_integral += matat_d(data->r, i_ph, i_th) * matat_d(model->P_lm_th, i_th, PT(l, m)) * cos_mph * vecat_d(sinth, i_th);
                s_integral += matat_d(data->r, i_ph, i_th) * matat_d(model->P_lm_th, i_th, PT(l, m)) * sin_mph * vecat_d(sinth, i_th);

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
}

data_iso *load_data_iso(const char *filename, int t, int p) {

    data_iso *data = (data_iso *) malloc(sizeof(*data));
    if (!data) err(1, "Cannot allocate data points structure");

    printf("[load_data_iso] Opening %s\n", filename);

    const double d_th = PI / t; 
    const double d_ph = TWO_PI / p;

	data->N = t * p;
    data->p = p;
    data->t = t;
    data->dt = d_th;
    data->dp = d_ph;

    // Create the data object matrix
    data->th = Matrix_new_d(1, t);
    data->ph = Matrix_new_d(1, p);
    data->r  = Matrix_new_d(p, t); 

    // shorthand aliases:
    Matrix_d *th = data->th, *ph = data->ph, *r = data->r;


    // data->th = (double *) malloc(sizeof(*data->th) * t); 
    // data->ph = (double *) malloc(sizeof(*data->ph) * p);
	// data->r  = malloc(data->N * sizeof(double));

    printf("[load_data_iso] Trying to read %d elements\n", data->N);

	if (th == NULL || ph == NULL || r == NULL)
		err(1, "cannot allocate data points\n");
    

    // fill theta matrix
    // data->th[0] = 0;
    vecset_d(th, 0, 0);
    
    // vecset
    for (int i = 1; i < t; i++) {
        // data->th[i] = data->th->data[i - 1] + d_th;
        // vecset_d(th, i, vecat_d(th, i - 1) + d_th);
        *vecptr_d(th, i) = vecat_d(th, i - 1) + d_th;
        // data->data
    }

    // fill phi array 
    // data->ph[0] = PI;
    // vecset_d(ph, 0, PI);
    // for (int i = 1; i < p / 2; i++) {
    //     // data->ph[i] = data->ph->data[i - 1] + d_ph;
    //     *vecptr_d(ph, i) = vecat_d(ph, i - 1) + d_ph;
    // }
    // // data->ph[p / 2] = 0.0;
    // vecset_d(ph, p / 2, 0.0);
    // for (int i = (p / 2) + 1; i < p; i++) {
    //     // data->ph[i] = data->ph->data[i - 1] + d_ph;
    //     *vecptr_d(ph, i) = vecat_d(ph, i - 1) + d_ph;
    // }

    vecset_d(ph, 0, - PI);
    for (int i = 1; i < p; i++) {
        // data->ph[i] = data->ph->data[i - 1] + d_ph;
        *vecptr_d(ph, i) = vecat_d(ph, i - 1) + d_ph;
    }


	FILE *f = fopen(filename, "r");
	if (f == NULL)
		err(1, "cannot open %s", filename);

    double p_lambda;
    double p_phi;
    double p_val;

    int count = 0;
	for (int i = 0; i < data->N; i++) {

		// int k = fscanf(f, "%lg %lg %lg", &data->lambda[i], &data->phi[i], &data->V[i]);
		int k = fscanf(f, "%lg %lg %lg", &p_lambda, &p_phi, &p_val);

        // data->r[i] = p_val;
        vecset_d(r, i, p_val);
        count++;

		if (k == EOF) {
			if (ferror(f))
				err(1, "read error");
			errx(1, "premature end-of-file after %d records", i);
		}
		if (k != 3)
			errx(1, "parse error on line %d", i+1);
	}
	fclose(f);

    printf("[load_data_iso] Read %d lines\n", count);

    return data;

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
    double *C = cs_mat->data;
    // double *S = matacc_d(cs_mat, 1, 0); // halfway point

    int index = 0;;
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

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

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

// Compute the gradient
double compute_gradient_clm(const data_iso *data, const Matrix_d *clm, const Matrix_d *pcs, int l, int m) {

    double sum = 0;

    for (int i = 0; i < data->N; i++) {

        double clmi = matat_d(clm, i, PT(l, m));
        sum += 2 * clmi * pcs->data[i] - 2 * clmi * vecat_d(data->r, i);
    }

    return sum;
}

double compute_gradient_clm_points(const iso_model *iso, int l, int m, const Matrix_i *indices) {

    double sum = 0;

    // iterate through the indices
    for (size_t i = 0; i < Matrix_size_i(indices); i++) {
        int index = indices->data[i];
        double clmi = matat_d(iso->model->clm, index, PT(l, m));
        sum += 2 * clmi * iso->model->pcs->data[index] - 2 * clmi * vecat_d(iso->data->r, index);
    }

    return sum;
}

double compute_gradient_slm_points(const iso_model *iso, int l, int m, const Matrix_i *indices) {

    double sum = 0;

    // iterate through the indices
    for (size_t i = 0; i < Matrix_size_i(indices); i++) {
        int index = indices->data[i];
        double slmi = matat_d(iso->model->slm, index, PT(l, m));
        sum += 2 * slmi * iso->model->pcs->data[index] - 2 * slmi * vecat_d(iso->data->r, index);
    }

    return sum;
}

double compute_gradient_slm(const data_iso *data, const Matrix_d *slm, const Matrix_d *pcs, int l, int m) {


    double sum = 0;

    for (int i = 0; i < data->N; i++) {

        double slmi = matat_d(slm, i, PT(l, m));
        // sum += 2 * slmi * pcs->data[i] - 2 * slmi * data->r[i];
        sum += 2 * slmi * pcs->data[i] - 2 * slmi * vecat_d(data->r, i);
    }

    return sum;
}

// Compute gradient in the form 
// [dMSE/dc00 dMSE/dc10 ... dMSE/dclm
//  dMSE/ds00 dMSE/ds10 ... dMSE/dslm]
// Assume model has C_lm and S_lm and P_lm_th
Matrix_d *compute_gradient(const data_iso *data, const spherical_model *model) {

    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    const Matrix_d *clm = NULL;
    const Matrix_d *slm = NULL;
    const Matrix_d *pcs = NULL;

    printf("Hi im in gradient\n");

    printf("%p, %p, %p\n", model->clm, model->slm, model->pcs);

    // check if the model has the models already computed
    if (model->clm != NULL && model->slm != NULL && model->pcs != NULL) {
    // We will actually need to recompute the indices for a general purpose calculation
    // if (false && model->clm != NULL && model->slm != NULL && model->pcs != NULL) {
        printf("Already computed coefficients\n");
        clm = model->clm; // 2 x ll matrix
        slm  = model->slm;
        pcs = model->pcs; // sum(clm * Plmcos * cos)_i for i in 1..N
    } else {
        clm = compute_mse_coeff_clm(data, model); // 1 x ll matrix
        slm  = compute_mse_coeff_slm(data, model); // 1 x ll matrix
        pcs = compute_pcs(data, model, clm, slm); // sum(clm * Plmcos * cos)_i for i in 1..N
    }

    // allocate space for output array

    // double *grad = (double *) malloc(sizeof(*grad) * 2 * ll);
    Matrix_d *grad = Matrix_new_d(2, ll);

    for (size_t l = 0; l <= model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {

            grad->data[PT(l, m)] = compute_gradient_clm(data, clm, pcs, l, m);
            grad->data[PT(l, m) + ll] = compute_gradient_slm(data, slm, pcs, l, m);

        }
    }

    return grad;
}

Matrix_d *compute_stochastic_gradient(const iso_model *iso, int n) {


    const int ll = iso->model->ll;

    const Matrix_d *clm = clm = iso->model->clm; // 2 x ll matrix
    const Matrix_d *slm = iso->model->slm;
    const Matrix_d *pcs = iso->model->pcs; 

    // generate random indices
    Matrix_i *indices = runif_i(n, 0, iso->data->N);

    printf("Hi im in gradient\n");

    // printf("%p, %p, %p\n", iso->model->clm, iso->model->slm, iso->model->pcs);
    Matrix_d *grad = Matrix_new_d(2, ll);

    for (size_t l = 0; l <= iso->model->lmax; l++) {
        for (size_t m = 0; m <= l; m++) {

            grad->data[PT(l, m)]      = compute_gradient_clm_points(iso, l, m, indices);
            grad->data[PT(l, m) + ll] = compute_gradient_slm_points(iso, l, m, indices);

        }
    }

    Matrix_free_i(indices);

    return grad;

}


// Assume that a model object has been fully initialized
double compute_mse(const iso_model *iso) {

    const spherical_model *model = iso->model;
    const data_iso        *data  = iso->data;


    const int ll = model->ll;

    double sum = 0;

    // Compute MSE by getting a list of predictions
    // const Matrix_d *predictions = modelPredictN(iso, iso->data->N);
    // const Matrix_d *predictions = modelPredictData(iso);
    const Matrix_d *predictions = modelPredictData(iso);

    // now compute the difference between the two
    for (int i = 0; i < data->N; i++) {
        double diff = vecat_d(predictions, i) - vecat_d(data->r, i);
        sum += diff * diff;
    }

    return sum / data->N;
}

double estimate_mse(const iso_model *iso, int n) {

    // estimate the MSE using only i data points
    // Matrix_i *indices = runif_i(n, 0, iso->data->N);
    Matrix_i *indices = Matrix_ij_i(1, n);

    // Matrix_print_i(indices);

    const spherical_model *model = iso->model;
    const data_iso        *data  = iso->data;

    double sum = 0;

    const Matrix_d *predictions = modelPredictDataPoints(iso, indices);

    for (int i = 0; i < n; i++) {
        int index = indices->data[i] - 1;
        double diff = vecat_d(predictions, i) - vecat_d(data->r, index);
        sum += diff * diff;
    }

    return sum / n;
}

double compute_average_error(const iso_model *iso) {

    const spherical_model *model = iso->model;
    const data_iso        *data  = iso->data;


    const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    double sum = 0;

    // Compute MSE by getting a list of predictions
    const Matrix_d *predictions = modelPredictN(iso, iso->data->N);

    // now compute the difference between the two
    for (int i = 0; i < data->N; i++) {
        sum += vecat_d(predictions, i) - vecat_d(data->r, i);
    }

    return sum / data->N;
}

// Given a gradient matrix and a model, tweak the values of 
// C_lm to improve the score
// gradient has size (N * 2) x ((lmax + 1)(lmax + 2)/2)
void adjust_parameters(const Matrix_d *grad, spherical_model *model, double alpha) {

    // const int ll = (model->lmax + 1) * (model->lmax + 2) / 2;

    // create a scaled version of the gradient vector (not necessary, but since the gradient vector should 
    // be reasonably small, we can perform this operation with very little cost)
    Matrix_d *grad_scaled = Matrix_clone_d(grad); // FAST memory copy
    matmultscalar_d(grad_scaled, -alpha); // Just multiply the first row
    // Matrix_mult_k

    // loop through the gradient matrix and adjust the model
    // printf("Gradient: \n");
    // Matrix_print_d(grad);


    // printf("Negative gradient * alpha\n");
    // Matrix_print_d(grad_scaled);

    // matadd_d(model->)
    MatIter_row_add_row_d(Matrix_row_begin_d(model->C_lm, 0), 
                          Matrix_row_end_d(model->C_lm, 0), 
                          Matrix_row_begin_d(grad_scaled, 0));

    MatIter_row_add_row_d(Matrix_row_begin_d(model->S_lm, 0),
                          Matrix_row_end_d(model->S_lm, 0),
                          Matrix_row_begin_d(grad_scaled, 1));
   
    // for (int l = 0; l <= model->lmax; l++) {
    //     for (int m = 0; m <= l; m++) {

    //         // printf ("C(%d, %d): %lf\n", l, m, model->C_lm[PT(l, m)]);
    //        model_C_plus(model, l, m, -alpha * matat_d(grad, 0, PT(l, m))) -alpha * matat_d(grad, 0, PT(l, m));
    //        model_S(model, l, m) += -alpha * matat_d(grad, 1, PT(l, m));
    //     }
    // }
}