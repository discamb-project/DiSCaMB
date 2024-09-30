SET(TARGET_H 
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/algebra3d.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/CrossProductLcs.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/geometry3d.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/GraphAlgorithms.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/InjectionsIterator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/LebedevGrid.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/MathUtilities.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/Matrix3.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/MixedRadixNumberIterator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/molecular_radial_grids.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/NaturalCubicSpline.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/PeriodicGridInterpolator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/RadialGrid.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/RealSphericalHarmonics.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/set_theory.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/SphConverter.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/SphMatrix.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/statistics.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/Vector2.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/Vector3.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/SphericalHarmonics.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/Polynomial.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/MathUtilities/Polynomial.h" )
    
#SET(TARGET_SRC
#    "${CMAKE_CURRENT_SOURCE_DIR}/geometry3d.cpp")      

SET(TARGET_SRC
    algebra3d.cpp
    CrossProductLcs.cpp 
    geometry3d.cpp
	GraphAlgorithms.cpp
	InjectionsIterator.cpp
    LebedevGrid.cpp
    MixedRadixNumberIterator.cpp
    molecular_radial_grids.cpp
    NaturalCubicSpline.cpp 
    PeriodicGridInterpolator.cpp
    RadialGrid.cpp
    RealSphericalHarmonics.cpp 
    SphConverter.cpp
    SphericalHarmonics.cpp
    SphMatrix.cpp
    statistics.cpp
    )


SETUP_OBJ_LIB(MathUtilities)