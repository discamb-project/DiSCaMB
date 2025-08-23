IF(DISCAMB_BUILD_DOCS)
    find_package(Doxygen REQUIRED)

    if(DOXYGEN_FOUND)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
        add_custom_target(doc 
                          ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                          COMMENT "Generating API documentation with Doxygen" VERBATIM
                             )
        
        set_target_properties(doc 
            PROPERTIES 
                EXCLUDE_FROM_ALL TRUE
                FOLDER "doc")
        
        MESSAGE("DOC PROCESSED")
    else()
        MESSAGE("Doxygen not found, API documentations will not be generated")
    endif()

ENDIF(DISCAMB_BUILD_DOCS)
