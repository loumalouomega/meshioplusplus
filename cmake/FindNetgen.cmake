# FindNetgen.cmake
# ----------------
# Locate nglib, the library interface of the Netgen mesher
# (https://github.com/NGSolve/netgen, LGPL-2.1), the backend of the
# tetrahedralize operation (roadmap 6.7). meshio++ calls nglib's C-style
# functions only (Ng_Init, Ng_AddPoint, Ng_AddSurfaceElement,
# Ng_GenerateVolumeMesh, ...) and links the shared library: Netgen is
# bring-your-own (never vendored, never auto-downloaded), as KaHIP is.
#
# Search order
# 1. NETGEN_ROOT / NETGEN_DIR environment variables (or CMake cache variables)
# 2. Standard system paths (Debian and Ubuntu ship libnglib-dev)
#
# Result variables
# ----------------
# NETGEN_FOUND        - TRUE if nglib is found
# NETGEN_INCLUDE_DIRS - The directory holding nglib.h (and its include/)
# NETGEN_LIBRARIES    - libnglib and libngcore
#
# Targets
# -------
# netgen::nglib       - IMPORTED INTERFACE target linking both

set(_netgen_search_hints)
foreach(_v NETGEN_ROOT NETGEN_DIR)
    if(DEFINED ENV{${_v}})
        list(APPEND _netgen_search_hints "$ENV{${_v}}")
    endif()
    if(DEFINED ${_v})
        list(APPEND _netgen_search_hints "${${_v}}")
    endif()
endforeach()

# Debian puts nglib.h in include/netgen/, a source install in include/.
find_path(NETGEN_INCLUDE_DIR
    NAMES nglib.h
    HINTS ${_netgen_search_hints}
    PATH_SUFFIXES include/netgen netgen include
    DOC "Directory containing nglib.h")
find_library(NETGEN_NGLIB_LIBRARY
    NAMES nglib
    HINTS ${_netgen_search_hints}
    PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu
    DOC "The nglib library")
find_library(NETGEN_NGCORE_LIBRARY
    NAMES ngcore
    HINTS ${_netgen_search_hints}
    PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu
    DOC "Netgen's core library, which nglib needs")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Netgen
    REQUIRED_VARS NETGEN_NGLIB_LIBRARY NETGEN_NGCORE_LIBRARY NETGEN_INCLUDE_DIR)

if(NETGEN_FOUND)
    set(NETGEN_INCLUDE_DIRS "${NETGEN_INCLUDE_DIR}")
    if(EXISTS "${NETGEN_INCLUDE_DIR}/include")
        list(APPEND NETGEN_INCLUDE_DIRS "${NETGEN_INCLUDE_DIR}/include")
    endif()
    set(NETGEN_LIBRARIES "${NETGEN_NGLIB_LIBRARY};${NETGEN_NGCORE_LIBRARY}")
    if(NOT TARGET netgen::nglib)
        add_library(netgen::nglib INTERFACE IMPORTED)
        set_target_properties(netgen::nglib PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${NETGEN_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${NETGEN_LIBRARIES}")
    endif()
endif()

mark_as_advanced(NETGEN_INCLUDE_DIR NETGEN_NGLIB_LIBRARY NETGEN_NGCORE_LIBRARY)
