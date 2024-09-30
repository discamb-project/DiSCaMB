include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../third-party/ucif")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../third-party/ucif/antlr3")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../third-party/ucif/antlr3/include")

SET(TARGET_H "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/atom_type_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/cif_io.h" 
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/discamb_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/fragmentation_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/gaussian_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/hkl_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/MATTS_BankReader.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/mol2_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/molecule_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/olex2_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/NativeIAM_Reader.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/pdb_io.h" 
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/proatom_db_io.h" 
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/rsp_io.h" 
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/shelx_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/structure_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/tham_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/tsc_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/vtk_io.h" 
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/wfn_io.h" 
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/xd_io.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/IO/xyz_io.h")
             
SET(TARGET_SRC 
    atom_type_io.cpp
    cif_io.cpp
    discamb_io.cpp
    fragmentation_io.cpp
    gaussian_io.cpp
    hkl_io.cpp
    MATTS_BankReader.cpp
    mol2_io.cpp
    molecule_io.cpp
    olex2_io.cpp
    NativeIAM_Reader.cpp
    pdb_io.cpp
    proatom_db_io.cpp
    rsp_io.cpp
    shelx_io.cpp
    structure_io.cpp
    tham_io.cpp
    tsc_io.cpp
    vtk_io.cpp
    wfn_io.cpp
    xd_io.cpp
    xyz_io.cpp)

#include_directories(${HDF5_INCLUDE_DIRS})             

#list(APPEND TARGET_H "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb /IO/EdIntegrationGridIO.h")
#list(APPEND TARGET_SRC "EdIntegrationGridIO.cpp")

ADD_LIBRARY(IO OBJECT ${TARGET_H} ${TARGET_SRC})

SET_PROPERTY(TARGET IO PROPERTY FOLDER "libs")
if(MT_MSVC_RUNTIME_LIB AND MSVC)
    set_property(TARGET IO PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif(MT_MSVC_RUNTIME_LIB AND MSVC) 

