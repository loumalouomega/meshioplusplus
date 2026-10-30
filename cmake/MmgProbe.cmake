# What the found MMG was built with, read from its library and headers (MMG's own
# config records neither). Included by the root build and by cmake/addon_mmg.cmake
# once Mmg::libmmg_so exists. Sets:
#   _mio_mmg_lib        the libmmg shared library file
#   _mio_mmg_incs       its include directories
#   _mio_mmg_int64      non-empty when MMG5_int is 64-bit
#   _mio_mmg_elas       non-empty when libmmg references ELAS (LS_elastic)
#   _mio_mmg_libdir     the directory of libmmg
#   _mio_mmg_elas_libs  libElas and libCommons by path (libmmg records no RPATH
#                       to them), empty without ELAS
#   _mio_mmg_version    mmg_VERSION, or "unknown"
# The library file and the headers, read for what MMG was built with: the
# width of MMG5_int (roadmap 3.2.5.2) from mmgcmakedefines.h, and ELAS -- the
# linear-elasticity library move needs -- from libmmg's reference to its
# LS_elastic, since MMG's headers do not record USE_ELAS.
set(_mio_mmg_lib "")
get_target_property(_mio_mmg_configs Mmg::libmmg_so IMPORTED_CONFIGURATIONS)
foreach(_c IN LISTS _mio_mmg_configs)
  get_target_property(_mio_mmg_lib Mmg::libmmg_so IMPORTED_LOCATION_${_c})
  if(_mio_mmg_lib)
    break()
  endif()
endforeach()
get_target_property(_mio_mmg_incs Mmg::libmmg_so INTERFACE_INCLUDE_DIRECTORIES)
set(_mio_mmg_int64 "")
foreach(_inc IN LISTS _mio_mmg_incs)
  if(EXISTS "${_inc}/mmg/common/mmgcmakedefines.h")
    file(STRINGS "${_inc}/mmg/common/mmgcmakedefines.h" _mio_mmg_int64
         REGEX "define MMG5_int int64_t")
    break()
  endif()
endforeach()
set(_mio_mmg_elas "")
if(_mio_mmg_lib AND EXISTS "${_mio_mmg_lib}")
  file(STRINGS "${_mio_mmg_lib}" _mio_mmg_elas REGEX "LS_elastic" LIMIT_COUNT 1)
endif()

set(_mio_mmg_elas_libs "")
if(_mio_mmg_lib)
get_filename_component(_mio_mmg_libdir "${_mio_mmg_lib}" DIRECTORY)
endif()
if(_mio_mmg_elas)
find_library(MESHIOPLUSPLUS_ELAS_LIBRARY NAMES Elas HINTS "${_mio_mmg_libdir}" ENV ELAS_DIR
             PATH_SUFFIXES lib)
find_library(MESHIOPLUSPLUS_ELAS_COMMONS_LIBRARY NAMES Commons
             HINTS "${_mio_mmg_libdir}" ENV ELAS_DIR PATH_SUFFIXES lib)
foreach(_l MESHIOPLUSPLUS_ELAS_LIBRARY MESHIOPLUSPLUS_ELAS_COMMONS_LIBRARY)
  if(${_l})
    list(APPEND _mio_mmg_elas_libs "${${_l}}")
  endif()
endforeach()
endif()
set(_mio_mmg_version "${mmg_VERSION}")
if(NOT _mio_mmg_version)
set(_mio_mmg_version "unknown")
endif()
