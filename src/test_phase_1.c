#include "geodesy.h"

int main() {

    const int lmax = 5;
    const int n_theta = 20;
    const double d_theta = PI / n_theta;
    const char datafile[] = "small";

    // Let's create a theta matrix 
    Matrix_d *theta = Matrix_new_d(1, n_theta);

    for (int i = 0; i < n_theta; i++) {
        theta->data[i] = d_theta * i;
    }

    Matrix_print_d(theta);

    write_binary_plm(lmax, theta, datafile);

    char binary_filename[50] = {0};
    sprintf(binary_filename, "ETOPO1_%s_P%d.bin", datafile, lmax);

    printf("Trying to read from %s\n", binary_filename);

    Matrix_d *plm = read_binary_plm(lmax, n_theta, binary_filename);

    Matrix_print_d(plm);

    // Damn man, this binary file exchange looks PRETTY FUCKING GOOD

    Matrix_d *p0 = read_binary_plm_l(lmax, 0, n_theta, binary_filename);
    Matrix_d *p1 = read_binary_plm_l(lmax, 1, n_theta, binary_filename);
    Matrix_d *p2 = read_binary_plm_l(lmax, 2, n_theta, binary_filename);
    Matrix_d *p3 = read_binary_plm_l(lmax, 3, n_theta, binary_filename);

    Matrix_print_d(p0);
    Matrix_print_d(p1);
    Matrix_print_d(p2);
    Matrix_print_d(p3);

    // Again, just fucking working LET's GO

    return 0;
}