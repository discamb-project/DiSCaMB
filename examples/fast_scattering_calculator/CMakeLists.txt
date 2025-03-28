add_example_exec(fast_scattering_calculator 
               fast_scattering_calculator.cpp AnyScattererStructureFactorCalculator2.cpp AnyScattererStructureFactorCalculator2.h)

#target_link_libraries(discamb_ff ${MPI_CXX_LIBRARIES})

target_link_libraries(fast_scattering_calculator discamb OpenMP::OpenMP_CXX)
