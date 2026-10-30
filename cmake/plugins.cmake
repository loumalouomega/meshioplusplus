# Provider plugins for the native CLI and the C API (roadmap 4.7.3), built with
# -DMESHIOPLUSPLUS_BUILD_PLUGINS=ON: one shared library per provider that is found,
# exporting the two entry points providers.hpp describes. A plugin is built
# WITHOUT linking the core: the host (the CLI, linked -rdynamic, or the C API
# library, promoted to global scope when it loads one) supplies the core's
# functions, so the plugin and the host share one registry. It needs the host's
# compiler, C++ runtime and release, which the ABI-version check makes explicit.
#
#   MESHIOPLUSPLUS_PLUGIN_PATH=<dir>:<dir>/meshioplusplus_plugin_occt.so meshioplusplus to-cad ...
#
# Windows is not supported: a DLL resolves its imports at link time against a named
# import library, so a plugin cannot borrow its host's symbols (the Python add-ons
# cover Windows).
if(WIN32)
  message(FATAL_ERROR "meshio++: MESHIOPLUSPLUS_BUILD_PLUGINS is not supported on Windows")
endif()

set(MESHIOPLUSPLUS_PLUGIN_DIR "${CMAKE_BINARY_DIR}/plugins" CACHE PATH
    "Where the plugin libraries are built")

# Defined once per provider: a MODULE from the provider's sources and its plugin/ entry points.
function(meshioplusplus_add_plugin name)
  file(GLOB _sources src/cpp/providers/${name}/*.cpp src/cpp/providers/${name}/plugin/*.cpp)
  add_library(meshioplusplus_plugin_${name} MODULE ${_sources})
  set_target_properties(meshioplusplus_plugin_${name} PROPERTIES
    PREFIX "" OUTPUT_NAME meshioplusplus_plugin_${name} CXX_STANDARD 20
    CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON
    LIBRARY_OUTPUT_DIRECTORY "${MESHIOPLUSPLUS_PLUGIN_DIR}")
  # The mesh backend and parallel macros the host was built with decide the layout
  # of the types a provider passes across: the plugin takes the core's own.
  target_compile_definitions(meshioplusplus_plugin_${name} PRIVATE
    MESHIOPLUSPLUS_MESH_BACKEND_${_meshioplusplus_mesh_backend} ${_mio_core_defs})
  target_include_directories(meshioplusplus_plugin_${name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/include ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/src)
  if(APPLE)
    # The core's functions are the host's, resolved at load.
    target_link_options(meshioplusplus_plugin_${name} PRIVATE -undefined dynamic_lookup)
  endif()
  install(TARGETS meshioplusplus_plugin_${name}
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/meshioplusplus/plugins)
endfunction()

# -- OpenCASCADE ---------------------------------------------------------------
find_package(OpenCASCADE CONFIG QUIET)
if(OpenCASCADE_FOUND AND NOT OpenCASCADE_VERSION VERSION_LESS 7.6)
  include(${CMAKE_CURRENT_LIST_DIR}/OcctToolkits.cmake)
  meshioplusplus_add_plugin(occt)
  target_compile_definitions(meshioplusplus_plugin_occt PRIVATE MESHIOPLUSPLUS_HAS_OCCT)
  target_include_directories(meshioplusplus_plugin_occt SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
  target_link_libraries(meshioplusplus_plugin_occt PRIVATE ${_mio_occt_libs})
  message(STATUS "meshio++: plugin meshioplusplus_plugin_occt (OpenCASCADE ${OpenCASCADE_VERSION})")
else()
  message(STATUS "meshio++: no OpenCASCADE 7.6 or later: no OpenCASCADE plugin")
endif()

# -- Netgen --------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
find_package(Netgen QUIET)
if(NETGEN_FOUND)
  meshioplusplus_add_plugin(netgen)
  target_compile_definitions(meshioplusplus_plugin_netgen PRIVATE MESHIOPLUSPLUS_HAS_NETGEN
    MESHIOPLUSPLUS_NETGEN_VERSION="${MESHIOPLUSPLUS_NETGEN_VERSION}")
  target_link_libraries(meshioplusplus_plugin_netgen PRIVATE netgen::nglib)
  message(STATUS "meshio++: plugin meshioplusplus_plugin_netgen")
else()
  message(STATUS "meshio++: no Netgen: no Netgen plugin")
endif()

# -- MMG (dynamic only) --------------------------------------------------------
find_package(mmg CONFIG QUIET)
if(TARGET Mmg::libmmg_so)
  include(${CMAKE_CURRENT_LIST_DIR}/MmgProbe.cmake)
  meshioplusplus_add_plugin(mmg)
  target_compile_definitions(meshioplusplus_plugin_mmg PRIVATE MESHIOPLUSPLUS_HAS_MMG
    MESHIOPLUSPLUS_MMG_VERSION="${_mio_mmg_version}")
  if(_mio_mmg_int64)
    target_compile_definitions(meshioplusplus_plugin_mmg PRIVATE MESHIOPLUSPLUS_MMG_INT64)
  endif()
  if(_mio_mmg_elas)
    target_compile_definitions(meshioplusplus_plugin_mmg PRIVATE MESHIOPLUSPLUS_HAS_MMG_ELAS)
  endif()
  target_link_libraries(meshioplusplus_plugin_mmg PRIVATE Mmg::libmmg_so ${_mio_mmg_elas_libs})
  message(STATUS "meshio++: plugin meshioplusplus_plugin_mmg")
else()
  message(STATUS "meshio++: no shared MMG: no MMG plugin")
endif()
