#include "geodesy.h"
#include "omp.h"

// Let's see how many threads are available on this machine

int main() {

    #pragma omp parallel
    {

        if (omp_get_thread_num() == 0) {
            printf("Num available threads: %d\n", omp_get_num_threads());
        }

    }

    return 0;
}