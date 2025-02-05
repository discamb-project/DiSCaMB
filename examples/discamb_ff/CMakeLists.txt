add_example_exec(discamb_ff
               discamb_ff.cpp )

#target_link_libraries(discamb_ff ${MPI_CXX_LIBRARIES})

target_link_libraries(discamb_ff discamb OpenMP::OpenMP_CXX)
