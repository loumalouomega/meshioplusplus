# The meshioplusplus-occt add-on's extension (MESHIOPLUSPLUS_ADDON=occt).
#
# Built by packages/python/meshioplusplus-occt against the installed base
# wheel: the core is found through the CMake package the wheel ships
# (meshioplusplus.get_cmake_dir()), never compiled again. The OpenCASCADE
# provider sources are the same ones a monolithic build compiles into the core.
# See doc/python_packages_plan.md.
include(${CMAKE_CURRENT_LIST_DIR}/addon_common.cmake)

find_package(OpenCASCADE CONFIG REQUIRED)
if(OpenCASCADE_VERSION VERSION_LESS 7.6)
  message(FATAL_ERROR "meshio++: OpenCASCADE ${OpenCASCADE_VERSION} is older than 7.6")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/OcctToolkits.cmake)

file(GLOB _mio_occt_sources src/cpp/providers/occt/*.cpp src/cpp/providers/occt/addon/*.cpp)
pybind11_add_module(_occt ${_mio_occt_sources})
target_compile_definitions(_occt PRIVATE MESHIOPLUSPLUS_HAS_OCCT)
target_include_directories(_occt PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/src)
target_include_directories(_occt SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries(_occt PRIVATE meshioplusplus::core_meshio ${_mio_occt_libs})
meshioplusplus_finish_addon(_occt meshioplusplus_occt)
