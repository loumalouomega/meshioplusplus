# The meshioplusplus-mmg add-on's extension (MESHIOPLUSPLUS_ADDON=mmg).
#
# Built by packages/python/meshioplusplus-mmg against the installed base wheel,
# like cmake/addon_occt.cmake. libmmg is linked dynamically only (LGPL-3.0-or-later),
# and the wheel is built against an MMG without ELAS (which would make the
# combination GPL-3.0), so `move` raises. The provider sources are the ones a
# monolithic build compiles into the core.
include(${CMAKE_CURRENT_LIST_DIR}/addon_common.cmake)

find_package(mmg CONFIG REQUIRED)
if(NOT TARGET Mmg::libmmg_so)
  message(FATAL_ERROR "meshio++: the add-on links MMG dynamically: build MMG with "
                      "BUILD_SHARED_LIBS=ON (build/build-mmg.sh) and set mmg_DIR")
endif()
include(${CMAKE_CURRENT_LIST_DIR}/MmgProbe.cmake)

file(GLOB _mio_mmg_sources src/cpp/providers/mmg/*.cpp src/cpp/providers/mmg/addon/*.cpp)
pybind11_add_module(_mmg ${_mio_mmg_sources})
target_compile_definitions(_mmg PRIVATE MESHIOPLUSPLUS_HAS_MMG
  MESHIOPLUSPLUS_MMG_VERSION="${_mio_mmg_version}")
if(_mio_mmg_int64)
  target_compile_definitions(_mmg PRIVATE MESHIOPLUSPLUS_MMG_INT64)
endif()
if(_mio_mmg_elas)
  target_compile_definitions(_mmg PRIVATE MESHIOPLUSPLUS_HAS_MMG_ELAS)
endif()
target_include_directories(_mmg PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/src)
target_link_libraries(_mmg PRIVATE meshioplusplus::core_meshio Mmg::libmmg_so ${_mio_mmg_elas_libs})
meshioplusplus_finish_addon(_mmg meshioplusplus_mmg)
