
SET(TARGET_H "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/basis_set_data.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/CrystalFragmentMultipoleCluster.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/CrystalFragmentWfnCalculationSettings.h"
			 "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/distributed_multipoles.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/ElectronDensityCalculator.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/ElectronDensityPartition.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/ElectronDensityPartitionType.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/fragmentation.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/GaussianRunner.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/HirshfeldPartition.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/MbisPartition.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/OrcaExecutionTimeEstimator.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/OrcaRunner.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/ProatomDB.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/PromoleculeElectronDensity.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/SphericalAtomicDensity.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/WaveFunctionCalculationData.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/QuantumChemistry/WaveFunctionFileFormat.h" )
             
SET(TARGET_SRC basis_set_data.cpp
               CrystalFragmentMultipoleCluster.cpp
			   distributed_multipoles.cpp
               ElectronDensityCalculator.cpp
               ElectronDensityPartition.cpp
               ElectronDensityPartitionType.cpp
               fragmentation.cpp
               GaussianRunner.cpp
               HirshfeldPartition.cpp
               MbisPartition.cpp
               OrcaExecutionTimeEstimator.cpp
               OrcaRunner.cpp
               ProatomDB.cpp
               PromoleculeElectronDensity.cpp
               SphericalAtomicDensity.cpp
			   WaveFunctionCalculationData.cpp
               WaveFunctionDataGeneratorRunner.cpp)
            
ADD_LIBRARY(QuantumChemistry OBJECT ${TARGET_H} ${TARGET_SRC})

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(QuantumChemistry PUBLIC OpenMP::OpenMP_CXX)
endif()

SET_PROPERTY(TARGET QuantumChemistry PROPERTY FOLDER "libs") 
if(MT_MSVC_RUNTIME_LIB AND MSVC)
    set_property(TARGET QuantumChemistry PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif(MT_MSVC_RUNTIME_LIB AND MSVC) 

