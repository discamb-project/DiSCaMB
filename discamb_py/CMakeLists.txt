set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

find_package(pybind11 CONFIG REQUIRED)

# Create an extension module
add_library(discamb_py MODULE discamb_py.cpp CrystalStructure.cpp CrystalStructure.h har_utilities.cpp har_utilities.h taam_utilities.cpp taam_utilities.h)
target_link_libraries(discamb_py PUBLIC pybind11::module discamb OpenMP::OpenMP_CXX)

