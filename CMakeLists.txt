cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(geodesy)

enable_testing(ON)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_compile_options(
    -Werror
    $<$<CONFIG:Debug>:-g>
)

if(CMAKE_COMPILER_IS_GNUCC) 
    message("Using gcc, adding -fopenmp")
endif()

include_directories(inc)
add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(test)