#include "geodesy.h"

int main() {

    // This program runs in like 5 minutes so I wouldn't suggest using it more than once. The data is stored under
    // ${project_dir}/data/small

    // Working with the small data set, create default models from L = 0 to L = 20
    const int L = 0;
    const int LMAX = 20;
    const char* data_in = "../../csv/ETOPO1_small.csv";
    const int npoint = 64800;
    Clock *clock = Clock_new();

    char filename[50] = {0};

    for (int l = L; l <= LMAX; l++) {

        // adjust file name
        sprintf(filename, "old_model_%d_180_360.txt", l);

        Clock_tic(clock);
        old_model_quiet(l, npoint, data_in, filename);
        Clock_toc(clock);

        printf("Old model L: %d in %lfs\n", l, elapsed_time(clock));
    }

    return 0;
}