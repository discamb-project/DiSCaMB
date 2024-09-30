
include_directories("${CMAKE_CURRENT_SOURCE_DIR}//vflib-2.0.6/include")

SET(TARGET_H 
    "./vflib-2.0.6/src/sortnodes.h"
    "./vflib-2.0.6/include/allocpool.h"
    "./vflib-2.0.6/include/argedit.h"
    "./vflib-2.0.6/include/argloader.h"
    "./vflib-2.0.6/include/argraph.h"
    "./vflib-2.0.6/include/dict.h"
    "./vflib-2.0.6/include/error.h"
    "./vflib-2.0.6/include/gene.h"
    "./vflib-2.0.6/include/gene_mesh.h"
    "./vflib-2.0.6/include/match.h"
    "./vflib-2.0.6/include/sd_state.h"
    "./vflib-2.0.6/include/state.h"
    "./vflib-2.0.6/include/ull_state.h"
    "./vflib-2.0.6/include/ull_sub_state.h"
    "./vflib-2.0.6/include/vf_mono_state.h"
    "./vflib-2.0.6/include/vf_state.h"
    "./vflib-2.0.6/include/vf_sub_state.h"
    "./vflib-2.0.6/include/vf2_mono_state.h"
    "./vflib-2.0.6/include/vf2_state.h"
    "./vflib-2.0.6/include/vf2_sub_state.h"
    "./vflib-2.0.6/include/xsubgraph.h" 
     )

SET(TARGET_SRC 
    "./vflib-2.0.6/src/argedit.cc"
    "./vflib-2.0.6/src/argloader.cc"
    "./vflib-2.0.6/src/argraph.cc"
    "./vflib-2.0.6/src/error.cc"
    "./vflib-2.0.6/src/gene.cc"
    "./vflib-2.0.6/src/gene_mesh.cc"
    "./vflib-2.0.6/src/match.cc"
    "./vflib-2.0.6/src/sd_state.cc"
    "./vflib-2.0.6/src/sortnodes.cc"
    "./vflib-2.0.6/src/ull_state.cc"
    "./vflib-2.0.6/src/ull_sub_state.cc"
    "./vflib-2.0.6/src/vf_mono_state.cc"
    "./vflib-2.0.6/src/vf_state.cc"
    "./vflib-2.0.6/src/vf_sub_state.cc"
    "./vflib-2.0.6/src/vf2_mono_state.cc"
    "./vflib-2.0.6/src/vf2_state.cc"
    "./vflib-2.0.6/src/vf2_sub_state.cc"
    "./vflib-2.0.6/src/xsubgraph.cc" 
 )

ADD_LIBRARY(vflib OBJECT ${TARGET_H} ${TARGET_SRC})
SET_PROPERTY(TARGET vflib PROPERTY FOLDER "third-party")

if(MT_MSVC_RUNTIME_LIB AND MSVC)
    set_property(TARGET vflib PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif(MT_MSVC_RUNTIME_LIB AND MSVC)


 