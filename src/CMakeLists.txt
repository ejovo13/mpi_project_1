add_executable(mejovo model_ejovo.c)
add_executable(model model.c)
add_executable(validate validate.c)
add_executable(tp1 test_phase_1.c)

target_link_libraries(mejovo geodesy)
target_link_libraries(model geodesy)
target_link_libraries(validate geodesy)
target_link_libraries(tp1 geodesy)