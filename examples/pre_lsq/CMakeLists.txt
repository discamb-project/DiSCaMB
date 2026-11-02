add_example_exec(pre_lsq pre_lsq.cpp)

if(OpenMP_FOUND)
	target_link_libraries(pre_lsq OpenMP::OpenMP_CXX)              
endif()

