cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(geodesy)

add_compile_options(
    -Werror
    $<$<CONFIG:Debug>:-g>
)

include_directories(inc)
add_subdirectory(lib)
add_subdirectory(src)