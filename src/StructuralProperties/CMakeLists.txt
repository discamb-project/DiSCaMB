include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../third-party/vflib/vflib-2.0.6/include")

SET(TARGET_H 
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/atom_selection.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/AtomInFragmentScoreCalculator.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/ConnectivityAlgorithm.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/CovalentRadiousBondDetector.h" 
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/GenericConnectivityAlgorithm.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/MacromolecularStructuralInformation.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/MolecularDisorder.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/RingCalculator.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/SimpleAIF_ScoreCalculator.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/StructuralProperties/structural_properties.h")

SET(TARGET_SRC 
    atom_selection.cpp
    AtomInFragmentScoreCalculator.cpp 
    ConnectivityAlgorithm.cpp
    CovalentRadiousBondDetector.cpp
    MacromolecularStructuralInformation.cpp
    MolecularDisorder.cpp
    RingCalculator.cpp
    SimpleAIF_ScoreCalculator.cpp
    structural_properties.cpp)
  

SETUP_OBJ_LIB(StructuralProperties)