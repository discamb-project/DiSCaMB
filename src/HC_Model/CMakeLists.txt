SET(TARGET_H 
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/HC_Model/ClementiRaimondiData.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/HC_Model/ClementiRoettiData.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/HC_Model/DeformationValenceParameters.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/HC_Model/HC_AtomTypeParameters.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/HC_Model/HC_ModelParameters.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/HC_Model/hc_model_utilities.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/HC_Model/HC_WfnData.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/HC_Model/XdLocalCoordinateSystem.h")

SET(TARGET_SRC 
    "ClementiRaimondiData.cpp"
    "ClementiRoettiData.cpp"
    "DeformationValenceParameters.cpp"
    "hc_model_utilities.cpp"
    "HC_WfnData.cpp"
    "XdLocalCoordinateSystem.cpp" )

SETUP_OBJ_LIB(HC_Model)

  