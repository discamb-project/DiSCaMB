add_example_exec(typeStat
               typeStat.cpp )

if(OpenMP_FOUND)
	target_link_libraries(typeStat discamb OpenMP::OpenMP_CXX)              
endif()

