function(add_geodesy_test TEST)
    set(SRC test_${TEST}.c)
    set(EXE test_${TEST})
    add_executable(${EXE} ${SRC})
    target_link_libraries(${EXE} geodesy)
    add_test(NAME ${TEST}
             COMMAND ${EXE})
endfunction()

add_geodesy_test(i_to_lm)
add_geodesy_test(omp_proc_available)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/run_proc.sh 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})