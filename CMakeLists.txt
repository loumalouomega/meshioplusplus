cmake_minimum_required(VERSION 3.15...3.30)

project(
  meshioplusplus_core
  LANGUAGES C CXX
  DESCRIPTION "C++ core for the meshio++ mesh I/O library")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# The pybind11 extension is the default consumer of meshioplusplus_core_obj,
# but non-Python consumers (the GoogleTest suite, the Emscripten/WASM target
# below) must be configurable without ever locating a Python interpreter or
# pybind11 -- e.g. under emcmake, where neither exists for the wasm32 target.
option(MESHIOPLUSPLUS_BUILD_PYTHON "Build the pybind11 _core extension" ON)
if(MESHIOPLUSPLUS_BUILD_PYTHON)
  # scikit-build-core provides the right Python; locate the module-only component.
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)
endif()

# ZLIB is optional too: without it the VTU zlib compression path falls back to
# Python (whose zlib is always available in the stdlib). This keeps Windows CI
# and wheels buildable with no system libraries.
option(MESHIOPLUSPLUS_WITH_ZLIB "Build the C++ VTU zlib compression path" ON)
if(MESHIOPLUSPLUS_WITH_ZLIB AND NOT EMSCRIPTEN)
  find_package(ZLIB QUIET)
endif()
# Under Emscripten there is no system zlib to find_package() -- the bundled
# port (-sUSE_ZLIB=1) supplies both <zlib.h> and the implementation, but only
# once that flag reaches the *compiler* invocation (not just the linker), so
# it's applied directly to meshioplusplus_core_obj below instead of via
# find_package/ZLIB::ZLIB.

# Optional heavy dependencies. When absent, the corresponding format sources
# compile to empty translation units (#ifdef-guarded) and the Python
# implementations serve as the runtime fallback.
option(MESHIOPLUSPLUS_WITH_HDF5 "Build the HDF5-backed formats (CGNS, HMF, H5M, MED, XDMF-HDF)" ON)
option(MESHIOPLUSPLUS_WITH_NETCDF "Build the netCDF-backed formats (Exodus)" ON)

if(MESHIOPLUSPLUS_WITH_HDF5)
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
if(MESHIOPLUSPLUS_WITH_NETCDF)
  # Also look in ~/.local for a user-built netcdf-c (no-sudo installs).
  find_package(netCDF CONFIG QUIET HINTS $ENV{HOME}/.local/lib/cmake/netCDF)
  if(NOT netCDF_FOUND)
    find_library(NETCDF_LIBRARY netcdf PATHS $ENV{HOME}/.local/lib)
    find_path(NETCDF_INCLUDE_DIR netcdf.h PATHS $ENV{HOME}/.local/include)
  endif()
endif()

# --------------------------------------------------------------------------
# meshioplusplus_core_obj: the pybind11-free C++ core (format readers/writers + pugixml).
# It is compiled once and shared by the Python extension (`_core`) and, when
# MESHIOPLUSPLUS_BUILD_TESTS is on, the standalone GoogleTest binary. All include
# dirs / feature flags / optional-library links are attached PUBLIC here so
# both consumers inherit them.
# --------------------------------------------------------------------------
file(GLOB_RECURSE MESHIOPLUSPLUS_CORE_SOURCES
    cpp/src/*.cpp
    cpp/third_party/pugixml/pugixml.cpp)

add_library(meshioplusplus_core_obj OBJECT ${MESHIOPLUSPLUS_CORE_SOURCES})
set_target_properties(meshioplusplus_core_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(meshioplusplus_core_obj PUBLIC cxx_std_20)
target_include_directories(
  meshioplusplus_core_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
                         ${CMAKE_CURRENT_SOURCE_DIR}/cpp/third_party/pugixml)

if(EMSCRIPTEN)
  # Emscripten disables JS-catchable C++ exceptions by default (an uncaught
  # exception calls abort() instead of unwinding into a JS Error) -- readers/
  # writers here throw ReadError/WriteError routinely (e.g. on an unsupported
  # construct), and bindings_js/js_bindings.cpp relies on those surfacing as
  # normal catchable JS exceptions. -fwasm-exceptions enables the native Wasm
  # exception-handling proposal (supported by all current browsers/Node),
  # which is faster than the legacy JS-longjmp emulation and needs no
  # separate -sDISABLE_EXCEPTION_CATCHING= linker flag. Applied PUBLIC so both
  # the compile step (here) and whatever links these objects inherit it.
  target_compile_options(meshioplusplus_core_obj PUBLIC "-fwasm-exceptions")
  target_link_options(meshioplusplus_core_obj PUBLIC "-fwasm-exceptions")
endif()

# --------------------------------------------------------------------------
# gcov/lcov instrumentation for the coverage CI job. OFF by default, so a
# normal `pip install` never pays for it. Applied PUBLIC so every consumer of
# these objects (_core, meshioplusplus_tests) is instrumented and emits into
# the same .gcno/.gcda set -- that is what lets one lcov capture cover the C++
# core as exercised by BOTH the GoogleTest binary and the pytest suite.
# -fprofile-update=atomic is required: meshioplusplus::parallel_for runs the
# instrumented loops on several threads, and the default non-atomic counter
# updates race and silently under-count.
# --------------------------------------------------------------------------
option(MESHIOPLUSPLUS_COVERAGE "Instrument the C++ core for gcov/lcov coverage" OFF)
if(MESHIOPLUSPLUS_COVERAGE)
  if(NOT (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang"))
    message(FATAL_ERROR "meshio++: MESHIOPLUSPLUS_COVERAGE needs GCC or Clang, got ${CMAKE_CXX_COMPILER_ID}")
  endif()
  message(STATUS "meshio++: coverage instrumentation ON (--coverage, -O0)")
  target_compile_options(meshioplusplus_core_obj PUBLIC
                         --coverage -O0 -g -fprofile-update=atomic)
  target_link_options(meshioplusplus_core_obj PUBLIC --coverage)
endif()

if(MESHIOPLUSPLUS_WITH_ZLIB AND EMSCRIPTEN)
  message(STATUS "meshio++: VTU zlib compression enabled (Emscripten -sUSE_ZLIB=1 port)")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_HAS_ZLIB)
  target_compile_options(meshioplusplus_core_obj PUBLIC "-sUSE_ZLIB=1")
  target_link_options(meshioplusplus_core_obj PUBLIC "-sUSE_ZLIB=1")
elseif(MESHIOPLUSPLUS_WITH_ZLIB AND ZLIB_FOUND)
  message(STATUS "meshio++: VTU zlib compression enabled")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_HAS_ZLIB)
  target_link_libraries(meshioplusplus_core_obj PUBLIC ZLIB::ZLIB)
else()
  message(STATUS "meshio++: zlib not found/disabled - VTU zlib falls back to Python")
endif()

if(MESHIOPLUSPLUS_WITH_HDF5 AND HDF5_FOUND)
  message(STATUS "meshio++: HDF5 formats enabled (HDF5 ${HDF5_VERSION})")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_HAS_HDF5)
  target_include_directories(meshioplusplus_core_obj PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(meshioplusplus_core_obj PUBLIC ${HDF5_C_LIBRARIES})
  if(HDF5_IS_PARALLEL AND MPI_C_FOUND)
    target_include_directories(meshioplusplus_core_obj PUBLIC ${MPI_C_INCLUDE_DIRS})
    target_link_libraries(meshioplusplus_core_obj PUBLIC MPI::MPI_C)
  endif()
else()
  message(STATUS "meshio++: HDF5 not found/disabled - HDF5 formats fall back to Python")
endif()

if(MESHIOPLUSPLUS_WITH_NETCDF AND (netCDF_FOUND OR NETCDF_LIBRARY))
  message(STATUS "meshio++: netCDF formats enabled")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_HAS_NETCDF)
  if(netCDF_FOUND)
    target_link_libraries(meshioplusplus_core_obj PUBLIC netCDF::netcdf)
  else()
    target_include_directories(meshioplusplus_core_obj PUBLIC ${NETCDF_INCLUDE_DIR})
    target_link_libraries(meshioplusplus_core_obj PUBLIC ${NETCDF_LIBRARY})
  endif()
else()
  message(STATUS "meshio++: netCDF not found/disabled - Exodus falls back to Python")
endif()

# Eigen (header-only, vendored as a git submodule at cpp/third_party/eigen):
# used for the MED Fortran<->C transpose. Optional -- when the submodule is not
# checked out (e.g. an sdist without submodules) the code falls back to the
# hand-written transpose loop. Run `git submodule update --init` to enable it.
option(MESHIOPLUSPLUS_WITH_EIGEN "Use vendored Eigen for the MED transpose" ON)
if(MESHIOPLUSPLUS_WITH_EIGEN AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cpp/third_party/eigen/Eigen/Dense")
  message(STATUS "meshio++: Eigen enabled (MED transpose)")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_HAS_EIGEN)
  target_include_directories(
    meshioplusplus_core_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/cpp/third_party/eigen)
else()
  message(STATUS "meshio++: Eigen not found/disabled - MED transpose uses the plain loop")
endif()

# --------------------------------------------------------------------------
# Parallel backend for meshioplusplus::parallel_for (cpp/include/meshioplusplus/parallel.hpp).
# Selected at configure time; exactly one MESHIOPLUSPLUS_PARALLEL_* macro is defined:
#   AUTO   - (default) OpenMP if available, else STL, else SEQ. Prefers OpenMP
#            because it is portable and needs no TBB; STL without TBB would
#            silently run sequentially.
#   SEQ    - sequential (always available)
#   STL    - C++17 std::execution::par. Built into MSVC's STL; on
#            libstdc++ it requires TBB - probed below, with an automatic
#            fallback to SEQ (warning) when unusable (e.g. Apple libc++).
#   OPENMP - #pragma omp parallel for (hard error if OpenMP is missing)
#   TBB    - tbb::parallel_for        (hard error if TBB is missing)
# Adding a new backend (Kokkos, ...) = a branch here + one #elif block in
# cpp/include/meshioplusplus/parallel.hpp.
# --------------------------------------------------------------------------
set(MESHIOPLUSPLUS_PARALLEL_BACKEND "AUTO" CACHE STRING
    "Parallel backend for meshioplusplus::parallel_for: AUTO, SEQ, STL, OPENMP or TBB")
set_property(CACHE MESHIOPLUSPLUS_PARALLEL_BACKEND PROPERTY STRINGS AUTO SEQ STL OPENMP TBB)
string(TOUPPER "${MESHIOPLUSPLUS_PARALLEL_BACKEND}" _meshioplusplus_parallel)

# AUTO: prefer OpenMP (portable: libgomp on manylinux, MSVC built-in, libomp on
# macOS; needs no TBB), then the STL(+TBB) path, else the sequential backend.
if(_meshioplusplus_parallel STREQUAL "AUTO")
  find_package(OpenMP QUIET COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    set(_meshioplusplus_parallel "OPENMP")
  else()
    set(_meshioplusplus_parallel "STL")
  endif()
endif()

if(_meshioplusplus_parallel STREQUAL "STL")
  if(MSVC)
    target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_PARALLEL_STL)
  else()
    find_package(TBB CONFIG QUIET)
    include(CheckCXXSourceCompiles)
    if(TBB_FOUND)
      set(CMAKE_REQUIRED_LIBRARIES TBB::tbb)
    endif()
    check_cxx_source_compiles(
      "
      #include <algorithm>
      #include <execution>
      #include <vector>
      int main() {
        std::vector<int> v(4);
        std::for_each(std::execution::par, v.begin(), v.end(), [](int& x) { x = 1; });
        return v[0] - 1;
      }"
      MESHIOPLUSPLUS_HAS_PSTL)
    unset(CMAKE_REQUIRED_LIBRARIES)
    if(MESHIOPLUSPLUS_HAS_PSTL)
      target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_PARALLEL_STL)
      if(TBB_FOUND)
        target_link_libraries(meshioplusplus_core_obj PUBLIC TBB::tbb)
      endif()
    else()
      message(WARNING
        "meshio++: the STL parallel backend is unusable on this toolchain "
        "(libstdc++ needs TBB installed; Apple libc++ has no parallel STL) - "
        "falling back to the sequential backend. Install TBB or configure with "
        "-DMESHIOPLUSPLUS_PARALLEL_BACKEND=OPENMP.")
      set(_meshioplusplus_parallel "SEQ")
      target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_PARALLEL_SEQ)
    endif()
  endif()
elseif(_meshioplusplus_parallel STREQUAL "OPENMP")
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_PARALLEL_OPENMP)
  target_link_libraries(meshioplusplus_core_obj PUBLIC OpenMP::OpenMP_CXX)
elseif(_meshioplusplus_parallel STREQUAL "TBB")
  find_package(TBB CONFIG REQUIRED)
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_PARALLEL_TBB)
  target_link_libraries(meshioplusplus_core_obj PUBLIC TBB::tbb)
elseif(_meshioplusplus_parallel STREQUAL "SEQ")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_PARALLEL_SEQ)
else()
  message(FATAL_ERROR
    "meshio++: unknown MESHIOPLUSPLUS_PARALLEL_BACKEND '${MESHIOPLUSPLUS_PARALLEL_BACKEND}' "
    "(use SEQ, STL, OPENMP or TBB)")
endif()
message(STATUS "meshio++: parallel backend: ${_meshioplusplus_parallel}")

if(MESHIOPLUSPLUS_BUILD_PYTHON)
  # The pybind11 extension: bindings + the shared core object library.
  file(GLOB MESHIOPLUSPLUS_BINDING_SOURCES bindings/*.cpp)
  pybind11_add_module(_core ${MESHIOPLUSPLUS_BINDING_SOURCES})
  target_include_directories(_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bindings)
  target_link_libraries(_core PRIVATE meshioplusplus_core_obj)

  # Off by default: only manylinux wheel builds (paired with a newer-than-image
  # GCC, e.g. gcc-toolset-13 on manylinux_2_28) need this, to keep the .so's
  # GLIBCXX/CXXABI/GLIBC symbol requirements within the target policy's floor.
  # LDFLAGS env-var seeding of CMAKE_MODULE_LINKER_FLAGS was tried first and
  # silently didn't apply to this MODULE target, hence targeting it directly.
  option(MESHIOPLUSPLUS_STATIC_RUNTIME "Statically link libgcc/libstdc++ into _core (portable manylinux wheels)" OFF)
  if(MESHIOPLUSPLUS_STATIC_RUNTIME AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_options(_core PRIVATE -static-libgcc -static-libstdc++)
  endif()

  # Land the extension next to the pure-Python package inside the wheel.
  install(TARGETS _core DESTINATION meshioplusplus)
endif()

# --------------------------------------------------------------------------
# WebAssembly build (Emscripten + embind), producing the @meshioplusplus/wasm
# npm package's native artifact. EMSCRIPTEN is set automatically by CMake when
# configured through emcmake; MESHIOPLUSPLUS_BUILD_WASM is the explicit,
# documented opt-in on top of that (see build/configure-wasm.sh). This is
# entirely independent of MESHIOPLUSPLUS_BUILD_PYTHON -- a wasm configure runs
# with -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF since no Python/pybind11 exists for
# the wasm32 target, and links the same meshioplusplus_core_obj the Python
# extension and GoogleTest suite already share.
# --------------------------------------------------------------------------
option(MESHIOPLUSPLUS_BUILD_WASM "Build the Emscripten/embind JS bindings" ${EMSCRIPTEN})
if(MESHIOPLUSPLUS_BUILD_WASM)
  if(NOT EMSCRIPTEN)
    message(FATAL_ERROR
      "meshio++: MESHIOPLUSPLUS_BUILD_WASM requires the Emscripten toolchain "
      "(configure with emcmake, e.g. via build/configure-wasm.sh).")
  endif()
  file(GLOB MESHIOPLUSPLUS_JS_BINDING_SOURCES bindings_js/*.cpp)
  add_executable(meshioplusplus_wasm ${MESHIOPLUSPLUS_JS_BINDING_SOURCES})
  target_include_directories(meshioplusplus_wasm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bindings_js)
  target_link_libraries(meshioplusplus_wasm PRIVATE meshioplusplus_core_obj)
  # -sUSE_ZLIB=1 is not repeated here: meshioplusplus_core_obj already applies
  # it PUBLIC (both compile and link) when MESHIOPLUSPLUS_WITH_ZLIB is on, and
  # this target inherits it transitively via target_link_libraries above.
  target_link_options(
    meshioplusplus_wasm PRIVATE
    "--bind"
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sFORCE_FILESYSTEM=1"
    "-sEXPORTED_RUNTIME_METHODS=['FS']")
  set_target_properties(meshioplusplus_wasm PROPERTIES OUTPUT_NAME "meshioplusplus_wasm"
                                                        SUFFIX ".mjs")
endif()

# --------------------------------------------------------------------------
# Optional standalone C++ tests (GoogleTest via CTest). OFF by default so that
# `pip install .` / cibuildwheel never build them. Enable with
# `-DMESHIO_BUILD_TESTS=ON` for a direct CMake configure.
# --------------------------------------------------------------------------
option(MESHIOPLUSPLUS_BUILD_TESTS "Build the C++ GoogleTest suite" OFF)
if(MESHIOPLUSPLUS_BUILD_TESTS)
  include(CTest)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
    URL_HASH SHA256=7b42b4d6ed48810c5362c265a17faebe90dc2373c885e5216439d37927f02926)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  file(GLOB MESHIOPLUSPLUS_TEST_SOURCES cpp/tests/*.cpp)
  add_executable(meshioplusplus_tests ${MESHIOPLUSPLUS_TEST_SOURCES})
  target_include_directories(meshioplusplus_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/tests)
  target_link_libraries(meshioplusplus_tests PRIVATE meshioplusplus_core_obj GTest::gtest_main)

  include(GoogleTest)
  gtest_discover_tests(meshioplusplus_tests)
endif()
