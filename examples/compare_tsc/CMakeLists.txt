add_example_exec(compare_tsc compare_tsc.cpp)
if(OpenMP_FOUND)
	target_link_libraries(compare_tsc OpenMP::OpenMP_CXX)              
endif()

