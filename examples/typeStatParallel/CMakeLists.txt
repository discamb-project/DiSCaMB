add_example_exec(typeStatParallel
               typeStatParallel.cpp )

if(OpenMP_FOUND)
	target_link_libraries(typeStatParallel discamb OpenMP::OpenMP_CXX)              
endif()

