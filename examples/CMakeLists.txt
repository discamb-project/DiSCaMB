include(${DISCAMB_MODULES_DIR}/add_example_exec.cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build/examples)

ADD_SUBDIRECTORY(compare_adps)
ADD_SUBDIRECTORY(compare_tsc)
ADD_SUBDIRECTORY(example_02_form_factors)
#ADD_SUBDIRECTORY(fast_scattering_calculator)
ADD_SUBDIRECTORY(harBankMaker)
ADD_SUBDIRECTORY(lab)
ADD_SUBDIRECTORY(discamb_ff)
ADD_SUBDIRECTORY(taam_symmetry)
ADD_SUBDIRECTORY(typeStat)
ADD_SUBDIRECTORY(typeStatParallel)
ADD_SUBDIRECTORY(typeSubtype)







