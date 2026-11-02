add_example_exec(sf_calculator sf_calculator.cpp)

if(OpenMP_FOUND)
	target_link_libraries(sf_calculator OpenMP::OpenMP_CXX)              
endif()


