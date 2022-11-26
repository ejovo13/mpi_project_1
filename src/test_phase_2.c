#include "geodesy.h"



// Assume that phase 1 has already been computed. Specifically, assume that the file
// ETOPO1_ultra_P788.bin has been precomputed, with an L_max of 788

// We are going to go ahead and focus on ETOPO1_small_P1000.bin for the moment

int main() {

    // const int n_theta = 10800;

    const int lmax_u = 788;
    const int t_u = 10800;
    const int p_u = t_u * 2;

    printf("Ultra {lmax = %d, nt = %d, nu = %d\n", lmax_u, t_u, p_u);

    const int bin_lmax_s = 1000;
    const int t_s = 180;
    const int p_s = t_s * 2;

    // So currently I have a binary data set stored in ETOPO1_small.bin,
    // I have precomputed values of PLM in ETOPO1_small_P1000.bin,

    // I need to load those values and then compute a model, then store the model values as binary

    // We are only computing the coefficients for a SINGLE L value using the write_binary_cslm function. 
    // Let's test it out

    const int selected_l = 100;
    // 
    const char *p_lm_th_binary = "ETOPO1_small_P1000.bin";
    const char *binary_dataset = "ETOPO1_small.bin";

    data_iso *data = load_data_binary(binary_dataset, t_s, p_s, true);

    write_binary_cslm(selected_l, bin_lmax_s, data, p_lm_th_binary);








    return 0;
}