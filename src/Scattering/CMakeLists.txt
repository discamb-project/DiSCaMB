  
SET(TARGET_H 
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/AnyHcCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/AnyIamCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/AnyScattererStructureFactorCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/AnyScattererStructureFactorCalculator2.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/AtomicDensityTransform.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/AtomicFormFactorCalculationsManager.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/AtomRepresentativeInfo.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/AtomTypeHC_Parameters.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/BankHcFunctions.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/CombinedStructureFactorCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/ConstFormFactorCalculationsManager.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/gar_utilities.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppens_SF_Engine.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppens_SF_Engine2.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppens_SF_Engine3.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppensStructureFactorCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppensStructureFactorCalculatorDev.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HcAtomBankStructureFactorCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HcFormFactorCalculationsManager.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HirshfeldAtomModelSettings.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/IamFormFactorCalculationsManager.h" 
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/IamSfCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/multipole_scattering.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/NGaussianFormFactor.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/NGaussianFormFactorsTable.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/NumericalSphericalAtomFormFactor.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/statistics.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/scattering_utilities.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/ScatteringData.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/SF_CalcDataTypes.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/SF_Engine_DataTypes.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/SfCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/SlaterTypeOrbitalScattering.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/StockholderAtomBankFormFactorCalculationManager.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/StockholderAtomBankSfCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/StockholderAtomFormFactorCalcManager.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/StockholderAtomSfCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/taam_utilities.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/TscFileBasedSfCalculator.h")
    
  
    
SET(TARGET_SRC
    AnyHcCalculator.cpp
    AnyIamCalculator.cpp
    AnyScattererStructureFactorCalculator.cpp
    AnyScattererStructureFactorCalculator2.cpp    
    AtomicDensityTransform.cpp
    AtomicFormFactorCalculationsManager.cpp
    BankHcFunctions.cpp
    CombinedStructureFactorCalculator.cpp
    ConstFormFactorCalculationsManager.cpp
    ElectronFromXrayFormFactorCalculationsManager.cpp
    gar_utilities.cpp
    HansenCoppens_SF_Engine.cpp 
    HansenCoppens_SF_Engine2.cpp
    HansenCoppens_SF_Engine3.cpp
    HansenCoppensStructureFactorCalculator.cpp 
    # HansenCoppensStructureFactorCalculatorDev.cpp 
    HcAtomBankStructureFactorCalculator.cpp
    HcFormFactorCalculationsManager.cpp
    HirshfeldAtomModelSettings.cpp
    IamFormFactorCalculationsManager.cpp
    IamSfCalculator.cpp
    multipole_scattering.cpp 
    NGaussianFormFactorsTable.cpp 
    NGaussianFormFactor.cpp 
    NumericalSphericalAtomFormFactor.cpp
    scattering_utilities.cpp 
    SF_CalcDataTypes.cpp 
    SfCalculator.cpp
    SlaterTypeOrbitalScattering.cpp 
    statistics.cpp
    StockholderAtomBankFormFactorCalculationManager.cpp
    StockholderAtomBankSfCalculator.cpp
    StockholderAtomFormFactorCalcManager.cpp
    StockholderAtomSfCalculator.cpp
    taam_utilities.cpp
    TscFileBasedSfCalculator.cpp)


            
ADD_LIBRARY(Scattering OBJECT ${TARGET_H} ${TARGET_SRC})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Scattering PUBLIC OpenMP::OpenMP_CXX)
endif()

SET_PROPERTY(TARGET Scattering PROPERTY FOLDER "libs") 
if(MT_MSVC_RUNTIME_LIB AND MSVC)
    set_property(TARGET Scattering PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif(MT_MSVC_RUNTIME_LIB AND MSVC) 

