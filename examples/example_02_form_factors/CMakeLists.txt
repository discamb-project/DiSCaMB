add_example_exec(example_02_form_factors
               example_02.cpp  )


if(OpenMP_FOUND)
	target_link_libraries(example_02_form_factors OpenMP::OpenMP_CXX)
endif()