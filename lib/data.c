/**========================================================================
 * ?                          data.c
 * @brief   : Functions/Structures relating to the loading and processing
 *            of raw data (binary or text) to be converted into models 
 * @details : 
 * @author  : Evan Voyles
 * @email   : ejovo13@yahoo.com
 * @date    : 2022-10-15
 *========================================================================**/
#include "data.h"
// #include "geodesy.h"

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
			errx(1, "premature end-of-file after %d records", i);
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
    // data->r  = Matrix_new_d(p, t);

    data->r = (int16_t *) malloc(sizeof(*data->r) * t * p); 

    // shorthand aliases:
    Matrix_d *th = data->th, *ph = data->ph;
    
    int16_t *r = data->r;


    // data->th = (double *) malloc(sizeof(*data->th) * t); 
    // data->ph = (double *) malloc(sizeof(*data->ph) * p);
	// data->r  = malloc(data->N * sizeof(double));

    printf("[load_data_iso] Trying to read %d elements\n", data->N);

	if (th == NULL || ph == NULL || r == NULL)
		err(1, "cannot allocate data points\n");
    

    // fill theta matrix
    // data->th[0] = 0;
    vecset_d(th, 0, 0);
    for (int i = 1; i < t; i++) {
        *vecptr_d(th, i) = vecat_d(th, i - 1) + d_th;
    }

    vecset_d(ph, 0, - PI);
    for (int i = 1; i < p; i++) {
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

        data->r[i] = (int16_t) p_val;
        // vecset_d(r, i, p_val);
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
        fprintf(out, "%.15lf\t%.15lf\t%lf\n", vecat_d(data->th, i % data->t), vecat_d(data->ph, i / data->t), (float) data->r[i]);
    }

    fclose(out);

}

void load_spherical_harmonics(const char *filename, int lmax, struct spherical_harmonics *self)
{
	FILE *g = fopen(filename, "r");
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
		if (k != 4)
			errx(1, "parse error");
	}
	fclose(g);
}

void free_data_iso(data_iso *data) {

    if (data == NULL) return;

    if (data->th != NULL) 
        Matrix_free_d(data->th);

    if (data->ph != NULL)
        Matrix_free_d(data->ph);

    if (data->r != NULL) 
        free(data->r);

    free(data);

    // if 

    // if (data->l_indices != NULL)
    //     free(data->l_indices);

    // if (data->m_indices != NULL)
    //     free(data->m_indices);
}

void head_data(const data_iso *data) {

    printf("data_iso @%p: { .N = %d, .t = %d, .p = %d, \n",
        data, data->N, data->t, data->p);

    printf("\t\t      .th(1:10) = ");
    Vector_print_head_d(data->th, 10);
    printf("\t\t      .ph(1:10) = ");
    Vector_print_head_d(data->ph, 10);
    printf("\t\t      .r(1:10)  = |");

    for (int i = 0; i < 9; i++) {
        printf("%d ", data->r[i]);
    }
    printf("%d |\n", data->r[9]);
    printf("\t\t      .dt = %lf, .dp = %lf\n", data->dt, data->dp);
    printf("\t\t    }\n");
    // Matrix_print_d(data->th);

}