add_executable(mejovo model_ejovo.c)
add_executable(model model.c)
add_executable(validate validate.c)
add_executable(tp1 test_phase_1.c)
add_executable(tp2 test_phase_2.c)
add_executable(phase_1 phase_1.c)
add_executable(reduce_dataset reduce_dataset_size.c)

target_link_libraries(mejovo geodesy)
target_link_libraries(model geodesy)
target_link_libraries(validate geodesy)
target_link_libraries(tp1 geodesy)
target_link_libraries(tp2 geodesy)
target_link_libraries(phase_1 geodesy)
target_link_libraries(reduce_dataset geodesy)