add_example_exec(lab lab.cpp)
target_link_libraries(lab OpenMP::OpenMP_CXX)              
if (TARGET Eigen3::Eigen)
	target_link_libraries(lab Eigen3::Eigen)
endif (TARGET Eigen3::Eigen)


