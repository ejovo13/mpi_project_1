cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(geodesy)

enable_testing(ON)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_compile_options(
    -Werror
    $<$<CONFIG:Release>:-O>
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Debug>:-Wall>
    $<$<CONFIG:Debug>:-Wno-unknown-pragmas>
    # $<$<CONFIG:Debug>:-fsanitize=address>
    # $<$<CONFIG:Debug>:-static-libasan>
    # $<$<CONFIG:Release>:-g>
    # $<$<CONFIG:Release>:-Wall>
    # $<$<CONFIG:Release>:-Wno-unknown-pragmas>
    # $<$<CONFIG:Release>:-fsanitize=address>
    # $<$<CONFIG:Release>:-static-libasan>
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Debug build")
endif()

if(CMAKE_COMPILER_IS_GNUCC) 
    message("Using gcc, adding -fopenmp")
endif()

include_directories(inc)
add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(test)