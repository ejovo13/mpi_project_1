#include "geodesy.h"

void write_binary_plm(int lmax, const Matrix_d *th, const char *dataset) {

    const int LL = (lmax + 1) * (lmax + 2) / 2;
    char filename[100] = {0};

    sprintf(filename, "ETOPO1_%s_P%d.bin", dataset, lmax);

    // Allocate space for a single row of the matrix
    Matrix_d *plm = Matrix_new_d(1, LL);

    // initialize coefficients a and b
    struct spherical_harmonics sph_model;
    setup_spherical_harmonics(lmax, &sph_model);

    // Open up the binary file to write
    FILE *bin = fopen(filename, "wb");

    int n = Matrix_size_d(th);
    int n_bytes_per_row = sizeof(double) * LL;

    Clock *clock = Clock_new();

    Clock_tic(clock);
        for (int i = 0; i < n; i++) {
            // Now compute the p values
            computeP(&sph_model, plm->data, cos(th->data[i]));
            fwrite(plm->data, n_bytes_per_row, 1, bin);
        }
    Clock_toc(clock);

    fclose(bin);
    printf("Wrote to binary file: %s in %lfs\n", filename, elapsed_time(clock));

    free(clock);
}

Matrix_d *read_binary_plm(int lmax, int n_th, const char *binary_filename) {

    const size_t LL = (lmax + 1) * (lmax + 2) / 2;

    // Allocate a Matrix to store the contents
    Matrix_d *P_lm_th = Matrix_new_d(n_th, LL);

    size_t n_bytes_per_row = sizeof(double) * LL;

    FILE *bin = fopen(binary_filename, "rb");

    // now fill up this matrix, row by row
    for (int i = 0; i < n_th; i++) {
        fread(matacc_d(P_lm_th, i, 0), n_bytes_per_row, 1, bin);
    }

    fclose(bin);
    printf("Read from binary file: %s\n", binary_filename);

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
    printf("Read in %d x %d P_L_th matrix from %s\n", n_theta, l, binary_filename);

    free(clock);

    return P_L_th;
}