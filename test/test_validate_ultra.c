// Load in the ultra data set, compute a prediction and save it

#include "geodesy.h"

int main(int argc, char ** argv) {

    int this_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

    data_iso *data = get_data_ultra(true);
    printf("Ultra data loaded\n");

    // Now let's precomp values up to 600
    Precomp *precomp = newPrecomp(0, 600, 600, data, "test_p600.bin");

    // Load the coefficients from a data set
    SphericalModel *model = buildSphericalModelMPI(get_data_med(true), 600, 600, "sph_med_600.bin", false, false, 0, world_size, this_rank);

    printf("Model loaded\n");

    // Now go ahead and compute the predictions, f_hat...
    Matrix_f *f_hat = compute_prediction_mpi(model, precomp, data, "ultra", world_size, this_rank);

    save_prediction(f_hat, "test_ultra_f_hat.bin");

    MPI_Finalize();

}