MACRO(add_example_exec)
    add_executable(${ARGN})
    target_link_libraries(${ARGV0} discamb)

    SET_PROPERTY(TARGET ${ARGV0} PROPERTY CXX_STANDARD 17)
    SET_PROPERTY(TARGET ${ARGV0} PROPERTY FOLDER "Examples")
    if(MT_MSVC_RUNTIME_LIB AND MSVC)
        set_property(TARGET ${ARGV0} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
        message(STATUS "use MT for " ${ARGV0})
    endif(MT_MSVC_RUNTIME_LIB AND MSVC)
    
ENDMACRO(add_example_exec)

