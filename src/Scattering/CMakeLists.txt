
SET(TARGET_H
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/agreement_factors.h"
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
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/deformation_valence_calc.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/disordered_structure_fragments.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/FragHarMacromol.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/gar_utilities.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppens_SF_Engine.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppens_SF_Engine2.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppens_SF_Engine3.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppens_SF_EngineSymm.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppensStructureFactorCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppensStructureFactorCalculatorDev.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HcAtomBankStructureFactorCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HcAtomBankStructureFactorCalculator2.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HcFormFactorCalculationsManager.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HirshfeldAtomModelSettings.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/IamFormFactorCalculationsManager.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/IamSfCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/multipole_scattering.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/NGaussianFormFactor.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/NGaussianFormFactorsTable.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/NumericalDensityIamFfCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/NumericalDensityIamSfCalculator.h"
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
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/TaamSfCalculator.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/TaamSfCalculatorMultiOrderedImpl.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/TscFileBasedSfCalculator.h")



SET(TARGET_SRC
    agreement_factors.cpp
    AnyHcCalculator.cpp
    AnyIamCalculator.cpp
    AnyScattererStructureFactorCalculator.cpp
    AnyScattererStructureFactorCalculator2.cpp
    AtomicDensityTransform.cpp
    AtomicFormFactorCalculationsManager.cpp
    BankHcFunctions.cpp
    CombinedStructureFactorCalculator.cpp
    ConstFormFactorCalculationsManager.cpp
    deformation_valence_calc.cpp
    disordered_structure_fragments.cpp
    ElectronFromXrayFormFactorCalculationsManager.cpp
	FragHarMacromol.cpp
    gar_utilities.cpp
    HansenCoppens_SF_Engine.cpp
    HansenCoppens_SF_Engine2.cpp
    HansenCoppens_SF_Engine3.cpp
    HansenCoppens_SF_EngineSymm.cpp
    HansenCoppensStructureFactorCalculator.cpp
    # HansenCoppensStructureFactorCalculatorDev.cpp
    HcAtomBankStructureFactorCalculator.cpp
    HcAtomBankStructureFactorCalculator2.cpp
    HcFormFactorCalculationsManager.cpp
    HirshfeldAtomModelSettings.cpp
    IamFormFactorCalculationsManager.cpp
    IamSfCalculator.cpp
    multipole_scattering.cpp
    NGaussianFormFactorsTable.cpp
    NGaussianFormFactor.cpp
    NumericalDensityIamFfCalculator.cpp
    NumericalDensityIamSfCalculator.cpp
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
    TaamSfCalculator.cpp
    TaamSfCalculatorMultiOrderedImpl.cpp
    TscFileBasedSfCalculator.cpp)

if(TAAM_WITH_GPU_IMP_26)
	LIST(APPEND TARGET_H "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/Scattering/HansenCoppens_SF_Engine4.h")
	LIST(APPEND TARGET_SRC HansenCoppens_SF_Engine4.cpp)
endif(TAAM_WITH_GPU_IMP_26)

option(SYCL "Enable using SYCL for calculations.")
option(SYCL_TARGETS "SYCL targets to compile for.")
option(SYCL_INCLUDE_DIR "Directory to where sycl/sycl.hpp resides. Optional if SYCL headers are present in system directories.")
option(SYCL_LIBRARY_DIR "Directory to where libsycl.so resides. Optional if cmake can find it on its own.")

ADD_LIBRARY(Scattering OBJECT ${TARGET_H} ${TARGET_SRC})

if(SYCL)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-fsycl HAS_SYCL)
  if(NOT HAS_SYCL)
    message(SEND_ERROR "SYCL is not supported by the chosen compiler")
  endif()

  if(SYCL_INCLUDE_DIR)
    set_source_files_properties(HansenCoppens_SF_Engine4.cpp PROPERTIES INCLUDE_DIRECTORIES ${SYCL_INCLUDE_DIR})
  endif()

  set_source_files_properties(HansenCoppens_SF_Engine4.cpp PROPERTIES COMPILE_DEFINITIONS SYCL)

  if(SYCL_TARGETS)
    set_source_files_properties(HansenCoppens_SF_Engine4.cpp PROPERTIES COMPILE_FLAGS "-fsycl -fsycl-targets=${SYCL_TARGETS}")
    target_link_options(Scattering PUBLIC "-fsycl" "-fsycl-targets=${SYCL_TARGETS}")
  else()
    set_source_files_properties(HansenCoppens_SF_Engine4.cpp PROPERTIES COMPILE_FLAGS -fsycl)
    target_link_options(Scattering PUBLIC "-fsycl")
  endif()

  find_library(SYCL_LIB NAMES sycl HINTS "${SYCL_LIBRARY_DIR}" REQUIRED)
  target_link_libraries(Scattering PRIVATE "${SYCL_LIB}")
endif()

if(NOT SYCL)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(Scattering PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

SET_PROPERTY(TARGET Scattering PROPERTY FOLDER "libs")
if(MT_MSVC_RUNTIME_LIB AND MSVC)
  set_property(TARGET Scattering PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif(MT_MSVC_RUNTIME_LIB AND MSVC)

