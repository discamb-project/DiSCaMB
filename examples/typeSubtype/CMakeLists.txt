add_example_exec(typeSubtype
               typeSubtype.cpp )

if(OpenMP_FOUND)
	target_link_libraries(typeSubtype OpenMP::OpenMP_CXX)              
endif()

