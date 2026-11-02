add_example_exec(fragHarMacromol fragHarMacromol.cpp)

if(OpenMP_FOUND)
	target_link_libraries(fragHarMacromol OpenMP::OpenMP_CXX)              
endif()


