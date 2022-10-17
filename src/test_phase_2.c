#include "geodesy.h"



// Assume that phase 1 has already been computed. Specifically, assume that the file
// ETOPO1_ultra_P788.bin has been precomputed, with an L_max of 788

int main() {

    const int lmax = 788;
    const int n_theta = 10800;

    // need to load the data from a text file. That's not going to happen yet, we first need 
    // to improve our loading capabilities.

    // will need to write a function that will significantly reduce the size of the ultra data set.
    // From 25gb to like 0.4gb
    // Then we will actually be able to load 0.4 GB into memory

    // on

    // for 

    // write_binary_cslm()
    // Examine the contents of our binary data sest.

    // binary_dataset_to_csv(64800, "ETOPO1_small.bin", "ETOPO1_small_mod.csv");
    binary_dataset_to_csv(583200, "ETOPO1_med.bin", "ETOPO1_med_mod.csv");


    return 0;
}