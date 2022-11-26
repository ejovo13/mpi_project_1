add_subdirectory(libejovo)

add_library(
    geodesy 
    data.c
    harmonics.c 
    model.c 
    learn.c 
    phase_1.c 
    phase_2.c 
    old_model.c 
    phase_3.c 
    mpi_util.c 
    ranges.c
    cli.c
)

target_link_libraries(
    geodesy
    # asan
    ejovo 
    m 
    OpenMP::OpenMP_C 
    MPI::MPI_C
)