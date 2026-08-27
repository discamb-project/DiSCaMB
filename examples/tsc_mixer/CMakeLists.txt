add_example_exec(tsc_mixer
               tsc_mixer.cpp )

target_link_libraries(tsc_mixer discamb OpenMP::OpenMP_CXX)
