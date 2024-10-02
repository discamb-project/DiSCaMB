SET(TARGET_H "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/constants.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/discamb_env.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/discamb_version.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/file_system_utilities.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/on_error.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/Exception.h" 
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/parse_cmd.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/Sheduler.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/Task.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/Timer.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/Tribool.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/string_utilities.h"
             "${CMAKE_CURRENT_SOURCE_DIR}/../../include/discamb/BasicUtilities/utilities.h")

SET(TARGET_SRC discamb_env.cpp discamb_version.cpp file_system_utilities.cpp parse_cmd.cpp Sheduler.cpp Task.cpp Timer.cpp utilities.cpp)

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()


SETUP_OBJ_LIB(BasicUtilities)

