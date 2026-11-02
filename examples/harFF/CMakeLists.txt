add_example_exec(harFF
               harFF.cpp )

if(OpenMP_FOUND)
	target_link_libraries(harFF OpenMP::OpenMP_CXX)
endif()

