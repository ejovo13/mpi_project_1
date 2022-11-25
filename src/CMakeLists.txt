function(add_geodesy EXE) 
    set(SRC ${EXE}.c)
    add_executable(${EXE} ${SRC})
    target_link_libraries(${EXE} geodesy)
endfunction()

# add_geodesy(model_ejovo)
# add_geodesy(old_model)
add_geodesy(validate)
add_geodesy(test_phase_1)
add_geodesy(test_phase_2)
add_geodesy(test_phase_3)
# add_geodesy(phase_1)
# add_geodesy(phase_2)
# add_geodesy(reduce_dataset_size)
add_geodesy(benchmark)
add_geodesy(mse_analysis)
add_geodesy(new_model)
# add_geodesy(new_model_omp)
add_geodesy(new_model_mpi)
# add_geodesy(new_model_omp2)
add_geodesy(parallel_speedup)

macro(move_file FILE) 
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${FILE}
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endmacro()

move_file(build_small_models.sh)
move_file(build_med_models.sh)
move_file(build_hi_models.sh)
move_file(build_old_small.sh)
move_file(validate_old_small.sh)
move_file(benchmark_mpi.sh)
move_file(run_mpi.sh)
move_file(world_ascii.txt)

# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/build_small_models.sh
#      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    
# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/build_med_models.sh
#      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/build_hi_models.sh
#      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/build_old_small.sh
#      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/validate_old_small.sh
#      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/benchmark_mpi.sh
#      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/run_mpi.sh
#      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/run_mpi.sh
#      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# add_executable(mejovo model_ejovo.c)
# add_executable(old_model old_model.c)
# add_executable(validate validate.c)
# add_executable(tp1 test_phase_1.c)
# add_executable(tp2 test_phase_2.c)
# add_executable(tp3 test_phase_3.c)
# add_executable(phase_1 phase_1.c)
# add_executable(phase_2 phase_2.c)
# add_executable(reduce_dataset reduce_dataset_size.c)
# add_executable(benchmark benchmark.c)
# add_executable(mse_analysis mse_analysis.c)

# target_link_libraries(mejovo geodesy)
# target_link_libraries(old_model geodesy)
# target_link_libraries(validate geodesy)
# target_link_libraries(tp1 geodesy)
# target_link_libraries(tp2 geodesy)
# target_link_libraries(tp3 geodesy)
# target_link_libraries(phase_1 geodesy)
# target_link_libraries(phase_2 geodesy)
# target_link_libraries(reduce_dataset geodesy)
# target_link_libraries(benchmark geodesy)
# target_link_libraries(mse_analysis geodesy)