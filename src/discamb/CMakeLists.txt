add_library(discamb
            $<TARGET_OBJECTS:ucif>
            $<TARGET_OBJECTS:vflib>
            $<TARGET_OBJECTS:AtomTyping>
            $<TARGET_OBJECTS:BasicChemistry>
            $<TARGET_OBJECTS:BasicUtilities>
            $<TARGET_OBJECTS:CrystalStructure>  
            $<TARGET_OBJECTS:HC_Model>
            $<TARGET_OBJECTS:IO>
            $<TARGET_OBJECTS:MathUtilities>                                      
            $<TARGET_OBJECTS:QuantumChemistry>
            $<TARGET_OBJECTS:Scattering> 
            $<TARGET_OBJECTS:StructuralProperties>
            $<TARGET_OBJECTS:Test>            
            discamb.cpp)

SET_PROPERTY(TARGET discamb PROPERTY FOLDER "libs")
if(MT_MSVC_RUNTIME_LIB AND MSVC)
    set_property(TARGET discamb PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif(MT_MSVC_RUNTIME_LIB AND MSVC)

#IF(GENERATE_INSTALL)
    install(TARGETS discamb
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
message(STATUS "binary dir " ${CMAKE_CURRENT_BINARY_DIR})
message(STATUS "binary dir " ${COMPILE_PDB_OUTPUT_DIRECTORY})
message(STATUS "binary dir " ${TARGET_FILE_DIR})
install(FILES "$<TARGET_FILE_DIR:discamb>/discambd.pdb" DESTINATION lib OPTIONAL)

#ENDIF(GENERATE_INSTALL)

