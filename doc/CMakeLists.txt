IF(BUILD_DOCS)
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
                EXCLUDE_FROM_ALL TRUE)
        
        MESSAGE("DOC PROCESSED")
    endif(DOXYGEN_FOUND)

ENDIF(BUILD_DOCS)
