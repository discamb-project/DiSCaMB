SET(TARGET_H "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicChemistry/ChemicalElement.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicChemistry/chemical_element_data.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicChemistry/basic_chemistry_utilities.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicChemistry/MoleculeData.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicChemistry/periodic_table.h")
SET(TARGET_SRC "ChemicalElement.cpp"
               "chemical_element_data.cpp"
               "basic_chemistry_utilities.cpp"
               "periodic_table.cpp")

SETUP_OBJ_LIB(BasicChemistry)