
SET(TARGET_H 
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/AtomInCrystalID.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/ConstantLocalCoordinateSystem.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/Crystal.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/CrystallographicRational.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/crystal_structure_utilities.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/CrystalVarianceCovarianceMatrix.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/ReciprocalLatticeUnitCell.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/SpaceGroup.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/SpaceGroupOperation.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/StructuralParametersConverter.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/UnitCellContent.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/CrystalStructure/UnitCell.h")
             
SET(TARGET_SRC "AtomInCrystalID.cpp"
               "ConstantLocalCoordinateSystem.cpp"
               "Crystal.cpp"
               "CrystallographicRational.cpp"
               "crystal_structure_utilities.cpp"
               "CrystalVarianceCovarianceMatrix.cpp"
               "LocalCoordinateSystemInCrystal.cpp"
               "ReciprocalLatticeUnitCell.cpp"
               "SpaceGroup.cpp"
               "SpaceGroupOperation.cpp"
               "StructuralParametersConverter.cpp"
               "UnitCell.cpp"
               "UnitCellContent.cpp"
)

SETUP_OBJ_LIB(CrystalStructure)