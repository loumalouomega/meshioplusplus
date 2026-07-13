cmake_minimum_required(VERSION 3.15...3.30)

project(
  meshio_core
  LANGUAGES C CXX
  DESCRIPTION "C++ core for the meshio mesh I/O library")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# scikit-build-core provides the right Python; locate the module-only component.
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# ZLIB is optional too: without it the VTU zlib compression path falls back to
# Python (whose zlib is always available in the stdlib). This keeps Windows CI
# and wheels buildable with no system libraries.
option(MESHIO_WITH_ZLIB "Build the C++ VTU zlib compression path" ON)
if(MESHIO_WITH_ZLIB)
  find_package(ZLIB QUIET)
endif()

# Optional heavy dependencies. When absent, the corresponding format sources
# compile to empty translation units (#ifdef-guarded) and the Python
# implementations serve as the runtime fallback.
option(MESHIO_WITH_HDF5 "Build the HDF5-backed formats (CGNS, HMF, H5M, MED, XDMF-HDF)" ON)
option(MESHIO_WITH_NETCDF "Build the netCDF-backed formats (Exodus)" ON)

if(MESHIO_WITH_HDF5)
  find_package(HDF5 QUIET COMPONENTS C)
  if(NOT HDF5_FOUND)
    # Some distros ship only an MPI-flavoured HDF5 (e.g. Debian's
    # libhdf5-openmpi-dev); FindHDF5 prefers serial by default, so retry.
    set(HDF5_PREFER_PARALLEL ON)
    find_package(HDF5 QUIET COMPONENTS C)
  endif()
  # A parallel HDF5 needs mpi.h even for serial use of the API.
  if(HDF5_FOUND AND HDF5_IS_PARALLEL)
    find_package(MPI QUIET COMPONENTS C)
    if(NOT MPI_C_FOUND)
      set(HDF5_FOUND FALSE)  # unusable without MPI headers -> Python fallback
    endif()
  endif()
endif()
if(MESHIO_WITH_NETCDF)
  # Also look in ~/.local for a user-built netcdf-c (no-sudo installs).
  find_package(netCDF CONFIG QUIET HINTS $ENV{HOME}/.local/lib/cmake/netCDF)
  if(NOT netCDF_FOUND)
    find_library(NETCDF_LIBRARY netcdf PATHS $ENV{HOME}/.local/lib)
    find_path(NETCDF_INCLUDE_DIR netcdf.h PATHS $ENV{HOME}/.local/include)
  endif()
endif()

# C++ core sources (grows as formats are ported).
file(GLOB_RECURSE MESHIO_CORE_SOURCES
    bindings/*.cpp
    cpp/src/*.cpp
    cpp/third_party/pugixml/pugixml.cpp)

pybind11_add_module(_core ${MESHIO_CORE_SOURCES})
target_include_directories(
  _core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
                ${CMAKE_CURRENT_SOURCE_DIR}/bindings
                ${CMAKE_CURRENT_SOURCE_DIR}/cpp/third_party/pugixml)
target_compile_features(_core PRIVATE cxx_std_17)

if(MESHIO_WITH_ZLIB AND ZLIB_FOUND)
  message(STATUS "meshio: VTU zlib compression enabled")
  target_compile_definitions(_core PRIVATE MESHIO_HAS_ZLIB)
  target_link_libraries(_core PRIVATE ZLIB::ZLIB)
else()
  message(STATUS "meshio: zlib not found/disabled - VTU zlib falls back to Python")
endif()

if(MESHIO_WITH_HDF5 AND HDF5_FOUND)
  message(STATUS "meshio: HDF5 formats enabled (HDF5 ${HDF5_VERSION})")
  target_compile_definitions(_core PRIVATE MESHIO_HAS_HDF5)
  target_include_directories(_core PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(_core PRIVATE ${HDF5_C_LIBRARIES})
  if(HDF5_IS_PARALLEL AND MPI_C_FOUND)
    target_include_directories(_core PRIVATE ${MPI_C_INCLUDE_DIRS})
    target_link_libraries(_core PRIVATE MPI::MPI_C)
  endif()
else()
  message(STATUS "meshio: HDF5 not found/disabled - HDF5 formats fall back to Python")
endif()

if(MESHIO_WITH_NETCDF AND (netCDF_FOUND OR NETCDF_LIBRARY))
  message(STATUS "meshio: netCDF formats enabled")
  target_compile_definitions(_core PRIVATE MESHIO_HAS_NETCDF)
  if(netCDF_FOUND)
    target_link_libraries(_core PRIVATE netCDF::netcdf)
  else()
    target_include_directories(_core PRIVATE ${NETCDF_INCLUDE_DIR})
    target_link_libraries(_core PRIVATE ${NETCDF_LIBRARY})
  endif()
else()
  message(STATUS "meshio: netCDF not found/disabled - Exodus falls back to Python")
endif()

# Land the extension next to the pure-Python package inside the wheel.
install(TARGETS _core DESTINATION meshio)
