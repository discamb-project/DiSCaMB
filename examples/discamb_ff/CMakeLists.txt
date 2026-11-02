add_example_exec(discamb_ff
               discamb_ff.cpp )


target_link_libraries(discamb_ff discamb)
if(OpenMP_FOUND)
	target_link_libraries(discamb_ff discamb OpenMP::OpenMP_CXX)
endif()

