# The meshioplusplus-netgen add-on's extension (MESHIOPLUSPLUS_ADDON=netgen).
#
# Built by packages/python/meshioplusplus-netgen against the installed base
# wheel, like cmake/addon_occt.cmake. The provider sources are the ones a
# monolithic build compiles into the core. nglib is linked dynamically (LGPL-2.1).
# Set MESHIOPLUSPLUS_NETGEN_VERSION to the version built against; nglib.h carries none.
include(${CMAKE_CURRENT_LIST_DIR}/addon_common.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
find_package(Netgen REQUIRED)
set(MESHIOPLUSPLUS_NETGEN_VERSION "unknown" CACHE STRING
    "The Netgen version reported by the tetrahedral mesher provider")

file(GLOB _mio_netgen_sources src/cpp/providers/netgen/*.cpp
     src/cpp/providers/netgen/addon/*.cpp)
pybind11_add_module(_netgen ${_mio_netgen_sources})
target_compile_definitions(_netgen PRIVATE MESHIOPLUSPLUS_HAS_NETGEN
  MESHIOPLUSPLUS_NETGEN_VERSION="${MESHIOPLUSPLUS_NETGEN_VERSION}")
target_link_libraries(_netgen PRIVATE meshioplusplus::core_meshio netgen::nglib)
meshioplusplus_finish_addon(_netgen meshioplusplus_netgen)
