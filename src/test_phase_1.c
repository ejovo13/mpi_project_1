#include "geodesy.h"

void t_small_read_write();
void t_ultra_read();
void t_comp_small(int l_model);

int main() {

    // t_small_read_write();
    // t_ultra_read();

    // t_comp_small(5);
    // t_comp_small(10);
    // t_comp_small(20);
    // t_comp_small(30);
    // t_comp_small(50);
    // t_comp_small(100);
    // t_comp_small(150);
    // t_comp_small(175);
    // t_comp_small(190);
    // t_comp_small(195);
    // t_comp_small(196);
    // t_comp_small(197);
    // t_comp_small(198);
    // t_comp_small(199);
    // t_comp_small(200);
    // t_comp_small(201);
    // t_comp_small(250);
    // t_comp_small(300);
    // t_comp_small(400);
    t_comp_small(500);

    return 0;
}

void t_small_read_write() {

    const int lmax = 5;
    const int n_theta = 20;
    const double d_theta = PI / n_theta;
    const char dataset_resolution[] = "small";

    // Let's create a theta matrix 
    Matrix_d *theta = Matrix_new_d(1, n_theta);

    for (int i = 0; i < n_theta; i++) {
        theta->data[i] = d_theta * i;
    }

    Matrix_print_d(theta);


    char binary_filename[50] = {0};
    sprintf(binary_filename, "TEST_ETOPO1_%s_P%d.bin", dataset_resolution, lmax);

    printf("Trying to read from %s\n", binary_filename);

    write_binary_plm(lmax, theta, binary_filename);
    Matrix_d *plm = read_binary_plm(lmax, n_theta, binary_filename, true);

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

}

void t_ultra_read() {

    /**========================================================================
     *!                            Ultra Test
     *========================================================================**/
    // test loading in the entire Ultra_P788 data set
    const char *ultra_dataset = "ETOPO1_ultra_P788.bin";
    const int lmax_u = 788;
    const int ntheta_u = 10800;

    Clock *clock = Clock_new();

    Clock_tic(clock); 
    // Matrix_d *plm_ultra = read_binary_plm(lmax_u, ntheta_u, ultra_dataset);
    Matrix_d *plm_ultra = read_binary_plm_l(lmax_u, 500, ntheta_u, ultra_dataset);
    Clock_toc(clock);

    printf("Successfully loaded plm_ultra %lu x %lu matrix in %lfs\n", 
        plm_ultra->nrows, plm_ultra->ncols, elapsed_time(clock));


    // Matrix_print_row_d(plm_ultra, );
}

void t_comp_small(int l_model) {

    data_iso *data_small = get_data_small(true);
    head_data(data_small);

    // data_iso *data_med = get_data_med();
    // head_data(data_med);

    // Now that we have a data file, go ahead and create a Precomp object

    // for now, keep L0 as 0 by default
    const int L0 = 0;
    const int LF = l_model; // only load in the first l_model values of L

    const int LMAX = 1000; // max of the binary data set
    const char *plm_bin = "ETOPO1_small_P1000.bin";

    // const int LMAX = 50;
    // const char *plm_bin = "ETOPO1_small_P50.bin";

    // Using these precomputed values, go ahead and compute the coefficients, storing them
    // in a brand new model
    SphericalModel *model = NULL;
    Precomp *precomp = NULL;
    char coeff_file_bin[50] = {0};
    sprintf(coeff_file_bin, "sph_%d_small.bin", LF);

    // if the file already exists, then load it. If not, precompute it again.
    /**========================================================================
     *!                        Precomp 
     *========================================================================**/ 

    precomp = newPrecomp(L0, LF, LMAX, data_small, plm_bin);

    FILE *test_open = fopen(coeff_file_bin, "rb");
    if (test_open == NULL) {

        printf("%s not found, computing Clm and Slm coefficients\n", coeff_file_bin);
        model = newSphericalModel(LF);    
        double time = modelComputeCSlmPrecomp(model, data_small, precomp);
        printf("Computed coefficients for L = %d in %lfs\n", LF, time);
        SphericalModelToBIN(model, "small");

    } else {

        // Model file coeff_file_bin already exists, so load from it
        fclose(test_open);
        model = loadSphericalModel(coeff_file_bin, LF);

    }

    // freePrecomp(precomp);

    // Now go ahead and print the first 10 coefficients to the screen
    Vector_print_head_d(model->C_lm, 10);
    Vector_print_head_d(model->S_lm, 10);

    // Output to a text file for visual inspection
    SphericalModelToTXT(model, "small");

    freeSphericalModel(model);

    // Write the model to bin
    // Free the contents in RAM
    // freeSphericalModel(model);
    // model = loadSphericalModel(coeff_file_bin, LF);

    // Now read the file into memory
    // And write it to text file

}