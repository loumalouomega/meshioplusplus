# FindKaHIP.cmake
# ----------------
# Locate the KaHIP partitioning library (https://github.com/KaHIP/KaHIP, MIT).
# Ported from the Kratos KaHIPApplication's Find module (KratosMultiphysics/
# Kratos PR #14453), minus its ParHIP section: meshio++ links only the serial
# kaffpa interface, never the MPI-based ParHIP, so no MPI dependency exists.
#
# KaHIP is bring-your-own (never vendored, never auto-downloaded — the same
# find_package-prefix policy as Metis/MMG in Kratos): install it yourself and
# point this module at it.
#
# Search order
# 1. KAHIP_ROOT / KAHIP_DIR environment variables (or CMake cache variables)
# 2. Standard system paths
# 3. pkg-config (kahip; shipped by KaHIP >= 3.24 installs)
#
# Result variables
# ----------------
# KAHIP_FOUND               - TRUE if KaHIP is found
# KAHIP_INCLUDE_DIRS        - Directory containing kaHIP_interface.h
# KAHIP_LIBRARIES           - The serial libkahip library
#
# Targets
# -------
# kahip::kahip              - IMPORTED shared/static serial library

cmake_policy(SET CMP0057 NEW) # IN_LIST operator

# ── Hints from environment and cache ─────────────────────────────────────────
set(_kahip_search_hints)

if(DEFINED ENV{KAHIP_ROOT})
    list(APPEND _kahip_search_hints "$ENV{KAHIP_ROOT}")
endif()
if(DEFINED KAHIP_ROOT)
    list(APPEND _kahip_search_hints "${KAHIP_ROOT}")
endif()
if(DEFINED ENV{KAHIP_DIR})
    list(APPEND _kahip_search_hints "$ENV{KAHIP_DIR}")
endif()
if(DEFINED KAHIP_DIR)
    list(APPEND _kahip_search_hints "${KAHIP_DIR}")
endif()

# ── Header ────────────────────────────────────────────────────────────────────
find_path(KAHIP_INCLUDE_DIR
    NAMES kaHIP_interface.h
    HINTS ${_kahip_search_hints}
    PATH_SUFFIXES interface
    DOC "KaHIP include directory (containing kaHIP_interface.h)")

# ── Serial library ────────────────────────────────────────────────────────────
find_library(KAHIP_LIBRARY
    NAMES kahip kahip_static
    HINTS ${_kahip_search_hints}
    PATH_SUFFIXES lib lib64 build
    DOC "KaHIP serial library (libkahip.so / libkahip_static.a)")

# ── pkg-config fallback ───────────────────────────────────────────────────────
if(NOT KAHIP_LIBRARY OR NOT KAHIP_INCLUDE_DIR)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(PC_KAHIP QUIET kahip)
        if(PC_KAHIP_FOUND)
            if(NOT KAHIP_INCLUDE_DIR)
                set(KAHIP_INCLUDE_DIR "${PC_KAHIP_INCLUDE_DIRS}")
            endif()
            if(NOT KAHIP_LIBRARY)
                set(KAHIP_LIBRARY "${PC_KAHIP_LIBRARIES}")
            endif()
        endif()
    endif()
endif()

# ── Handle REQUIRED / QUIET ────────────────────────────────────────────────────
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KaHIP
    REQUIRED_VARS KAHIP_LIBRARY KAHIP_INCLUDE_DIR
    FAIL_MESSAGE "Could NOT find KaHIP. Set KAHIP_ROOT to the KaHIP install/deploy directory.")

# ── Populate result variables ─────────────────────────────────────────────────
if(KAHIP_FOUND)
    set(KAHIP_INCLUDE_DIRS "${KAHIP_INCLUDE_DIR}")
    set(KAHIP_LIBRARIES    "${KAHIP_LIBRARY}")

    # ── IMPORTED target: kahip::kahip ────────────────────────────────────────
    if(NOT TARGET kahip::kahip)
        add_library(kahip::kahip UNKNOWN IMPORTED)
        set_target_properties(kahip::kahip PROPERTIES
            IMPORTED_LOCATION             "${KAHIP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${KAHIP_INCLUDE_DIR}")
    endif()
endif()

mark_as_advanced(
    KAHIP_INCLUDE_DIR
    KAHIP_LIBRARY)
