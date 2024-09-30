include_directories("${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/antlr3")

SET(TARGET_H 
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3baserecognizer.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3basetree.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3basetreeadaptor.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3bitset.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3collections.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3commontoken.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3commontree.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3commontreeadaptor.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3commontreenodestream.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3convertutf.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3cyclicdfa.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3debugeventlistener.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3defs.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3encodings.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3errors.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3exception.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3filestream.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3input.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3interfaces.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3intstream.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3lexer.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3memory.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3parser.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3parsetree.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3recognizersharedstate.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3rewritestreams.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3string.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3tokenstream.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/include/antlr3treeparser.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/antlr3config.h" 
     "${CMAKE_CURRENT_SOURCE_DIR}/builder.h" 
     "${CMAKE_CURRENT_SOURCE_DIR}/cifLexer.h" 
     "${CMAKE_CURRENT_SOURCE_DIR}/cifParser.h" 
     "${CMAKE_CURRENT_SOURCE_DIR}/parser.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/utils.h"
     )

SET(TARGET_SRC 
    "${CMAKE_CURRENT_SOURCE_DIR}/antlr3/src/antlr3baserecognizer.c"
    "./antlr3/src/antlr3basetree.c"
    "./antlr3/src/antlr3basetreeadaptor.c"
    "./antlr3/src/antlr3bitset.c"
    "./antlr3/src/antlr3collections.c"
    "./antlr3/src/antlr3commontoken.c"
    "./antlr3/src/antlr3commontree.c"
    "./antlr3/src/antlr3commontreeadaptor.c"
    "./antlr3/src/antlr3commontreenodestream.c"
    "./antlr3/src/antlr3convertutf.c"
    "./antlr3/src/antlr3cyclicdfa.c"
    "./antlr3/src/antlr3encodings.c"
    "./antlr3/src/antlr3exception.c"
    "./antlr3/src/antlr3filestream.c"
    "./antlr3/src/antlr3inputstream.c"
    "./antlr3/src/antlr3intstream.c"
    "./antlr3/src/antlr3lexer.c"
    "./antlr3/src/antlr3parser.c"
    "./antlr3/src/antlr3rewritestreams.c"
    "./antlr3/src/antlr3string.c"
    "./antlr3/src/antlr3tokenstream.c"
    "./antlr3/src/antlr3treeparser.c"
    "cifLexer.cpp"
    "cifParser.cpp" 
 )

ADD_LIBRARY(ucif OBJECT ${TARGET_H} ${TARGET_SRC})
SET_PROPERTY(TARGET ucif PROPERTY FOLDER "third-party")
if(MT_MSVC_RUNTIME_LIB AND MSVC)
    set_property(TARGET ucif PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif(MT_MSVC_RUNTIME_LIB AND MSVC)


 