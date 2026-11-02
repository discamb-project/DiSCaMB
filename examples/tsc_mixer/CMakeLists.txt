add_example_exec(tsc_mixer
               tsc_mixer.cpp )

if(OpenMP_FOUND)
	target_link_libraries(tsc_mixer OpenMP::OpenMP_CXX)
endif()
