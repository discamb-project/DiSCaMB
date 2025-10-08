#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 3.15)

cmake_policy(SET CMP0091 NEW)

PROJECT(discamb)

option(DISCAMB_BUILD_DOCS "Build the documentation" OFF)
option(DISCAMB_BUILD_EXAMPLES "Build the example programs" OFF)
option(DISCAMB_MSVC_MT "An options for MSVC to build multi-threaded statically-linked version" OFF)

if(DISCAMB_MSVC_MT)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON) 

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(DISCAMB_MODULES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(MSVC)
  add_definitions(/MP)
endif()



if((CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
    option(USE_STD_SPH_BESSEL "use std::sph_bessel" OFF) 
else()
    option(USE_STD_SPH_BESSEL "use std::sph_bessel" ON) 
endif()

if(USE_STD_SPH_BESSEL)
	set(HAS_SPH_BESSEL true)
endif()


find_package(OpenMP)
set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_MINSIZEREL_POSTFIX "s")
set(CMAKE_RELWITHDEBINFO_POSTFIX "rd")


#SET(DISCAMB_INCLUDE_PATH "" CACHE PATH "dir with discamb headers folder 'discamb'")
#include_directories(${DISCAMB_INCLUDE_PATH})

if(WIN32)
    SET(CMAKE_DETECTED_WIN32 TRUE)
endif(WIN32)

# add the binary tree to the search path for include files
# so that we will find config.h

include_directories("${PROJECT_BINARY_DIR}/build/include")
INCLUDE_DIRECTORIES(BEFORE include third-party)
ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/src")


IF(DISCAMB_BUILD_EXAMPLES)
    ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/examples")    
ENDIF(DISCAMB_BUILD_EXAMPLES)


IF(DISCAMB_BUILD_DOCS)
    ADD_SUBDIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/doc")
ENDIF(DISCAMB_BUILD_DOCS)

#------ INSTALLATION ------

    IF(INSTALLATION_PREFIX) 
        SET(CMAKE_INSTALL_PREFIX  ${INSTALLATION_PREFIX})
    ELSE(INSTALLATION_PREFIX) 
        IF(CMAKE_SYSTEM_NAME STREQUAL "Windows")
            SET(CMAKE_INSTALL_PREFIX  ${CMAKE_BINARY_DIR}/build/)
        ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Windows")                    
    ENDIF(INSTALLATION_PREFIX) 
    
    install(DIRECTORY include/discamb DESTINATION include)
    install(FILES "${PROJECT_BINARY_DIR}/build/include/discamb/config.h" DESTINATION include/discamb)

#--------

configure_file ("${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/build/include/discamb/config.h")

set (DISCAMB_BUILD_DIR "${PROJECT_BINARY_DIR}/build")


