add_example_exec(taam_symmetry taam_symmetry.cpp)

if(OpenMP_FOUND)
	target_link_libraries(taam_symmetry OpenMP::OpenMP_CXX)              
endif()

