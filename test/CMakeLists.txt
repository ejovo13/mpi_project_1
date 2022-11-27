function(add_geodesy_test TEST)
    set(SRC test_${TEST}.c)
    set(EXE test_${TEST})
    add_executable(${EXE} ${SRC})
    target_link_libraries(${EXE} geodesy)
    target_compile_options(${EXE} PUBLIC -g)
    add_test(NAME ${TEST}
             COMMAND ${EXE})
endfunction()

add_geodesy_test(i_to_lm)
add_geodesy_test(omp_proc_available)
add_geodesy_test(range)
add_geodesy_test(hello_mpi)
add_geodesy_test(validate_ultra)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/run_proc.sh 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_SOURCE_DIR}/bin/models/sph_med_500.bin
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_SOURCE_DIR}/bin/models/sph_small_100.bin
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_SOURCE_DIR}/bin/models/sph_small_200.bin
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_SOURCE_DIR}/bin/models/sph_small_300.bin
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_SOURCE_DIR}/bin/models/sph_med_600.bin
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


# file(COPY ${CMAKE_SOURCE_DIR}/bin/models/sph_small_300.bin
#     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})