include(${DISCAMB_MODULES_DIR}/add_example_exec.cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/examples)

ADD_SUBDIRECTORY(compare_adps)
ADD_SUBDIRECTORY(compare_tsc)
ADD_SUBDIRECTORY(fast_scattering_calculator)
ADD_SUBDIRECTORY(lab)
ADD_SUBDIRECTORY(discamb_ff)
ADD_SUBDIRECTORY(test)

