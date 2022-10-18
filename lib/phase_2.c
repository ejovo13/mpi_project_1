#include "geodesy.h"

// Compute a full set of coefficients CLm and SLm, then store them to two binary files 
void write_binary_cslm(int L, int lmax, int n_theta, const data_iso *data, const char *binary_file_in) {

    const int ll = (lmax + 1) * (lmax + 2) / 2;

    Matrix_d *C_L = Matrix_new_d(1, L);
    Matrix_d *S_L = Matrix_new_d(1, L);

    // n_theta x (L + 1) matrix
    Matrix_d *P_L_th = read_binary_plm_l(lmax, L, n_theta, binary_file_in);

    printf("[ write_binary_cslm ] computing coefficients for L = %d \n", L);

    // Initialize integral values
    double c_integral = 0;
    double s_integral = 0;
    int count = 0;

    for (int m = 0; m <= L; m++) {

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

            double ph_iph = vecat_d(data->ph, i_ph);
            double cos_mph = cos(m * ph_iph); // could use a recursive relationship to calclute this faster
            double sin_mph = sin(m * ph_iph);

            // use the midpoint formula
            c_integral += data->r[i] * matat_d(P_L_th, i_th, m) * cos_mph * vecat_d(sinth, i_th);
            s_integral += data->r[i] * matat_d(P_L_th, i_th, m) * sin_mph * vecat_d(sinth, i_th);

            count ++;
        }

        c_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);
        s_integral *= (data->dp * data->dt) / (2.0 * TWO_PI);

        vecset_d(C_L, m, c_integral);
        vecset_d(S_L, m, s_integral);

        Matrix_free_d(sinth);

    }

    // Now that C_L and S_L have been computed, go ahead and write them to a binary file
    char c_fileout[50] = {0};
    char s_fileout[50] = {0};
    sprintf(c_fileout, "C_%05d.bin", L);
    sprintf(s_fileout, "S_%05d.bin", L);

    FILE *c_coeff_out = fopen(c_fileout, "wb");
    FILE *s_coeff_out = fopen(s_fileout, "wb");

    const int n_bytes = sizeof(double) * (L + 1);

    fwrite(C_L->data, n_bytes, 0, c_coeff_out);
    fwrite(S_L->data, n_bytes, 0, s_coeff_out);

    fclose(c_coeff_out);
    fclose(s_coeff_out);

    printf("Wrote C_L data to %s\n", c_fileout);
    printf("Wrote S_L data to %s\n", s_fileout);

    Matrix_free_d(C_L);
    Matrix_free_d(S_L);

}

int16_t *read_binary_dataset(int n_points, const char *binary_file_in) {

    // allocate space for matrix
    int16_t *z = (int16_t *) malloc(sizeof(*z) * n_points);

    const size_t n_bytes = sizeof(z) * n_points;

    // Now open up the file and read through it
    FILE *bin = fopen(binary_file_in, "rb");
    fread(z, n_bytes, 1, bin);

    return z;
}

// Write binary file containing the altitude values to a csv file for visual inspection
void write_binary_z_to_csv(int n_points, const int16_t *z, const char *output_csv) {

    FILE *out_csv = fopen(output_csv, "w");

    for (int i = 0; i < n_points; i++) {
        fprintf(out_csv, "%d\n", z[i]);
    }

    fclose(out_csv);
}

void binary_dataset_to_csv(int n_points, const char *binary_in, const char *output_csv) {

    int16_t *z = read_binary_dataset(n_points, binary_in);
    write_binary_z_to_csv(n_points, z, output_csv);
    free(z);

    printf("Converted binary file %s to %s\n", binary_in, output_csv);

}

void reduce_csv_to_binary(int n_points, const char *csv_in, const char *binary_out) {

    // open the data set and scan through
    FILE *f = fopen(csv_in, "r");
	if (f == NULL)
		err(1, "cannot open %s", csv_in);

    float trash_lambda;
    float trash_ph;
    float z_float;
    int16_t z_int;    

    const size_t n_bytes = sizeof(z_int);

    FILE *bin = fopen(binary_out, "wb");

	for (int i = 0; i < n_points; i++) {

		int k = fscanf(f, "%lf %lf %f\n", &trash_lambda, &trash_ph, &z_float);

        z_int = (int16_t) z_float;

        // write to binary file
        fwrite(&z_int, n_bytes, 1, bin);

		if (k == EOF) {
			if (ferror(f))
				err(1, "read error");
			errx(1, "premature end-of-file after %d records", i);
		}
		if (k != 3)
			errx(1, "parse error on line %d", i+1);

	}


    fclose(bin);
    fclose(f);

    printf("Walked through dataset: %s\n", csv_in);

}

data_iso *load_data_binary(const char *binary_file_in, int t, int p) {

    data_iso *data = (data_iso *) malloc(sizeof(*data));
    if (!data) err(1, "Cannot allocate data points structure");

    printf("[load_data_binary] Opening %s\n", binary_file_in);

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
    data->r  = (int16_t *) malloc(sizeof(*data->r) * t * p); 

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

	FILE *f = fopen(binary_file_in, "rb");
	if (f == NULL)
		err(1, "cannot open %s", binary_file_in);

    double p_lambda;
    double p_phi;
    double p_val;

    int count = 0;
	for (int i = 0; i < data->N; i++) {

		// int k = fscanf(f, "%lg %lg %lg", &data->lambda[i], &data->phi[i], &data->V[i]);
		int k = fscanf(f, "%lg %lg %lg", &p_lambda, &p_phi, &p_val);

        // data->r[i] = p_val;
        // vecset_d(r, i, p_val);
        data->r[i] = (int16_t) p_val;
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