add_library(discamb discamb.cpp)

target_link_libraries(discamb PUBLIC
            ucif
            vflib
            AtomTyping
            BasicChemistry
            BasicUtilities
            CrystalStructure
            HC_Model
            IO
            MathUtilities
            QuantumChemistry
            Scattering
            StructuralProperties
            Test)

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

