cmake_minimum_required(VERSION 3.15...3.30)

# VERSION must track pyproject.toml's `version`, conanfile.py's `version`,
# packages/vcpkg/meshioplusplus/vcpkg.json's `version`, src/wasm/package.json's
# `version` and the "../wasm" entry in src/viewer/package-lock.json (bump all
# six together on a release) -- it feeds mio_version(), the shared-library
# VERSION properties, and the find_package/pkg-config metadata of the C API.
project(
  meshioplusplus_core
  VERSION 8.4.0
  LANGUAGES C CXX
  DESCRIPTION "C++ core for the meshio++ mesh I/O library")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Off by default: only portable release artifacts (manylinux wheel builds --
# paired with a newer-than-image GCC, e.g. gcc-toolset-13 on manylinux_2_28 --
# and the standalone CLI binary, see cli.yml) need this, to keep the binary's
# GLIBCXX/CXXABI symbol requirements within the target policy's floor (Linux)
# or avoid a vcredist dependency (MSVC static CRT). Applied below to whichever
# of _core / meshioplusplus_cli is actually being built.
option(MESHIOPLUSPLUS_STATIC_RUNTIME
  "Statically link the C++ runtime (libgcc/libstdc++ on GNU, static CRT on MSVC)" OFF)

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

# Additional VTK XML block-compression codecs. Both default OFF: zlib is the
# only codec a portable build can assume, and it stays the write default, so a
# pure build reads and writes exactly what it always did. Never probed under
# Emscripten -- there is no -sUSE_ZSTD/-sUSE_LZ4 port, and a host library must
# not be picked up for a wasm target.
option(MESHIOPLUSPLUS_WITH_ZSTD "Build the VTK XML zstd compression path" OFF)
option(MESHIOPLUSPLUS_WITH_LZ4 "Build the VTK XML lz4 compression path" OFF)
if(MESHIOPLUSPLUS_WITH_ZSTD AND NOT EMSCRIPTEN)
  find_package(zstd QUIET)
  if(NOT zstd_FOUND)
    find_path(ZSTD_INCLUDE_DIR zstd.h)
    find_library(ZSTD_LIBRARY NAMES zstd libzstd)
  endif()
endif()
if(MESHIOPLUSPLUS_WITH_LZ4 AND NOT EMSCRIPTEN)
  find_package(lz4 QUIET)
  if(NOT lz4_FOUND)
    find_path(LZ4_INCLUDE_DIR lz4.h)
    find_library(LZ4_LIBRARY NAMES lz4 liblz4)
  endif()
endif()

# KaHIP (MIT), the optional quality backend of the partition operation. OFF by
# default so the default build stays dependency-free; the SFC method is the
# always-available fallback. Bring-your-own install located via the in-repo
# cmake/FindKaHIP.cmake (KAHIP_ROOT prefix / pkg-config) -- never vendored,
# never auto-downloaded. Only the serial kaffpa interface is linked (no ParHIP,
# hence no MPI). Never probed under Emscripten -- no port exists, and a host
# library must not be picked up for a wasm target.
option(MESHIOPLUSPLUS_WITH_KAHIP "Build the KaHIP partition backend" OFF)
if(MESHIOPLUSPLUS_WITH_KAHIP AND NOT EMSCRIPTEN)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  find_package(KaHIP QUIET)
endif()

# Optional heavy dependencies. When absent, the corresponding format sources
# compile to empty translation units (#ifdef-guarded) and the Python
# implementations serve as the runtime fallback.
option(MESHIOPLUSPLUS_WITH_HDF5 "Build the HDF5-backed formats (CGNS, HMF, H5M, MED, XDMF-HDF)" ON)
option(MESHIOPLUSPLUS_WITH_NETCDF "Build the netCDF-backed formats (Exodus)" ON)

if(MESHIOPLUSPLUS_WITH_HDF5)
  if(EMSCRIPTEN)
    # There is only ever one HDF5 for this target: the static wasm32 build
    # produced by build/build-wasm-deps.sh, which build/configure-wasm.sh puts
    # on CMAKE_FIND_ROOT_PATH. It is serial by construction, so the
    # MPI-flavour retry below has nothing to find and the MPI probe nothing to
    # guard against -- both are skipped rather than left to mis-fire. Module
    # mode is kept (not the installed hdf5-config.cmake) so the
    # HDF5_INCLUDE_DIRS/HDF5_C_LIBRARIES plumbing further down stays identical
    # for every platform.
    set(HDF5_USE_STATIC_LIBRARIES ON)
    find_package(HDF5 QUIET COMPONENTS C)
  else()
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
endif()
if(MESHIOPLUSPLUS_WITH_NETCDF)
  # Also look in ~/.local for a user-built netcdf-c (no-sudo installs).
  find_package(netCDF CONFIG QUIET HINTS $ENV{HOME}/.local/lib/cmake/netCDF)
  if(NOT netCDF_FOUND)
    find_library(NETCDF_LIBRARY netcdf PATHS $ENV{HOME}/.local/lib)
    find_path(NETCDF_INCLUDE_DIR netcdf.h PATHS $ENV{HOME}/.local/include)
  endif()
endif()

# Polyscope (MIT), the native viewer behind the CLI's `view`/`screenshot`
# verbs. OFF by default and, unlike every other optional dependency here,
# attached to the **CLI target only** -- never to meshioplusplus_core_obj. The
# core, the C API, the Fortran module and the wasm build therefore cannot see
# it at all, which is what keeps "the core is dependency-free" literally true
# rather than true-by-default: Polyscope pulls in OpenGL, GLFW and X11.
#
# Vendored as a git submodule at src/cpp/third_party/polyscope (it is in no package
# manager and its own docs offer only add_subdirectory/FetchContent). It has
# nested submodules, so `git submodule update --init --recursive` is required;
# with a plain --init the EXISTS guard below simply leaves the viewer out.
option(MESHIOPLUSPLUS_WITH_POLYSCOPE
       "Build the CLI's native viewer (needs OpenGL/GLFW; see src/cpp/third_party/polyscope)"
       OFF)

# --------------------------------------------------------------------------
# meshioplusplus_core_obj: the pybind11-free C++ core (format readers/writers + pugixml).
# It is compiled once and shared by the Python extension (`_core`) and, when
# MESHIOPLUSPLUS_BUILD_TESTS is on, the standalone GoogleTest binary. All include
# dirs / feature flags / optional-library links are attached PUBLIC here so
# both consumers inherit them.
# --------------------------------------------------------------------------
file(GLOB_RECURSE MESHIOPLUSPLUS_CORE_SOURCES
    src/cpp/src/*.cpp
    src/cpp/third_party/pugixml/pugixml.cpp)

add_library(meshioplusplus_core_obj OBJECT ${MESHIOPLUSPLUS_CORE_SOURCES})
set_target_properties(meshioplusplus_core_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(meshioplusplus_core_obj PUBLIC cxx_std_20)

# Several operations (interpolate's barycentric weighted sum, merge/partition's
# weighted accumulations, ...) are pinned byte-identical against a pure-Python/
# numpy twin, which never fuses a multiply-add into one rounding step. x86_64
# has no FMA instruction unless a target flag (-mfma/-march=...) opts in, so
# this was silently never an issue on the Linux/Windows CI legs -- but AArch64
# (e.g. GitHub's arm64 macos-latest runners) has FMA in baseline hardware, and
# Clang's default -ffp-contract=on permits fusing `acc = acc + w * v` within a
# single expression, changing the rounding by up to 1 ulp relative to numpy's
# separately-rounded multiply and add (this is exactly what broke
# tests/python/test_interpolate.py::test_cpp_matches_python on macOS). Disable
# contraction outright so this class of expression rounds identically on every
# target; there is no MSVC equivalent needed since /fp:precise (the default)
# never contracts.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(meshioplusplus_core_obj PUBLIC -ffp-contract=off)
endif()

# On MSVC, /MT vs /MD is a *compile*-time flag baked into each .obj (unlike
# GNU's -static-libgcc/-static-libstdc++, which is link-only) -- every object
# library feeding a static-runtime final link must be compiled with the same
# setting, or the link fails with unresolved __imp_-prefixed CRT externals
# (nearbyint/isnan/ecvt_s/... -- the dynamic-CRT objects expect an import lib
# that a static-CRT link doesn't provide). Never applied when building the
# Python extension: pybind11 modules must always match the interpreter's own
# (dynamic) CRT, and Windows wheels never set MESHIOPLUSPLUS_STATIC_RUNTIME.
if(MESHIOPLUSPLUS_STATIC_RUNTIME AND MSVC AND NOT MESHIOPLUSPLUS_BUILD_PYTHON)
  set_property(TARGET meshioplusplus_core_obj PROPERTY
    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_include_directories(
  meshioplusplus_core_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/include
                         ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/third_party/pugixml)

if(EMSCRIPTEN)
  # Emscripten disables JS-catchable C++ exceptions by default (an uncaught
  # exception calls abort() instead of unwinding into a JS Error) -- readers/
  # writers here throw ReadError/WriteError routinely (e.g. on an unsupported
  # construct), and bindings/wasm/js_bindings.cpp relies on those surfacing as
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

if(MESHIOPLUSPLUS_WITH_ZSTD AND (zstd_FOUND OR ZSTD_LIBRARY))
  message(STATUS "meshio++: VTK XML zstd compression enabled")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_HAS_ZSTD)
  if(zstd_FOUND)
    target_link_libraries(meshioplusplus_core_obj PUBLIC zstd::libzstd_shared)
  else()
    target_include_directories(meshioplusplus_core_obj PUBLIC ${ZSTD_INCLUDE_DIR})
    target_link_libraries(meshioplusplus_core_obj PUBLIC ${ZSTD_LIBRARY})
  endif()
elseif(MESHIOPLUSPLUS_WITH_ZSTD)
  message(STATUS "meshio++: zstd requested but not found - zstd files fall back to Python")
endif()

if(MESHIOPLUSPLUS_WITH_LZ4 AND (lz4_FOUND OR LZ4_LIBRARY))
  message(STATUS "meshio++: VTK XML lz4 compression enabled")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_HAS_LZ4)
  if(lz4_FOUND)
    target_link_libraries(meshioplusplus_core_obj PUBLIC LZ4::lz4_shared)
  else()
    target_include_directories(meshioplusplus_core_obj PUBLIC ${LZ4_INCLUDE_DIR})
    target_link_libraries(meshioplusplus_core_obj PUBLIC ${LZ4_LIBRARY})
  endif()
elseif(MESHIOPLUSPLUS_WITH_LZ4)
  message(STATUS "meshio++: lz4 requested but not found - lz4 files fall back to Python")
endif()

if(MESHIOPLUSPLUS_WITH_KAHIP AND KAHIP_FOUND)
  message(STATUS "meshio++: KaHIP partition backend enabled (${KAHIP_LIBRARIES})")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_HAS_KAHIP)
  target_link_libraries(meshioplusplus_core_obj PUBLIC kahip::kahip)
elseif(MESHIOPLUSPLUS_WITH_KAHIP)
  message(STATUS "meshio++: KaHIP requested but not found (set KAHIP_ROOT) - "
                 "partition method 'kahip' will raise")
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

# Eigen (header-only, vendored as a git submodule at src/cpp/third_party/eigen):
# used for the MED Fortran<->C transpose. Optional -- when the submodule is not
# checked out (e.g. an sdist without submodules) the code falls back to the
# hand-written transpose loop. Run `git submodule update --init` to enable it.
option(MESHIOPLUSPLUS_WITH_EIGEN "Use vendored Eigen for the MED transpose" ON)
if(MESHIOPLUSPLUS_WITH_EIGEN AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/third_party/eigen/Eigen/Dense")
  message(STATUS "meshio++: Eigen enabled (MED transpose)")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_HAS_EIGEN)
  target_include_directories(
    meshioplusplus_core_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/third_party/eigen)
else()
  message(STATUS "meshio++: Eigen not found/disabled - MED transpose uses the plain loop")
endif()

# --------------------------------------------------------------------------
# Parallel backend for meshioplusplus::parallel_for (src/cpp/include/meshioplusplus/parallel.hpp).
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
# src/cpp/include/meshioplusplus/parallel.hpp.
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

# ---------------------------------------------------------------------------
# In-memory mesh backend (src/cpp/include/meshioplusplus/mesh.hpp). Exactly one
# MESHIOPLUSPLUS_MESH_BACKEND_* macro is defined; all format code is written
# against the uniform mesh API (mesh_api.hpp) so any backend compiles:
#   MESHIO - (default) the meshio-mirroring Mesh/CellBlock over dtype-erased
#            NDArrays. REQUIRED when MESHIOPLUSPLUS_BUILD_PYTHON=ON (the
#            zero-copy numpy boundary in bindings/python/np_conversions.hpp is
#            written against it).
#   NATIVE - canonical statically-typed storage (Float64 points, Int64
#            connectivity, CellType enum, CSR ragged blocks). The fastest
#            pure-C++ consumer surface; the WebAssembly build uses it.
#   KRATOS - a Kratos-Multiphysics-style ModelPart (Nodes/Elements/
#            Conditions/SubModelParts) behind the same API, for near-costless
#            exchange with Kratos (see kratos_bridge.hpp).
# ---------------------------------------------------------------------------
set(MESHIOPLUSPLUS_MESH_BACKEND "MESHIO" CACHE STRING
    "In-memory mesh backend: MESHIO, NATIVE or KRATOS")
set_property(CACHE MESHIOPLUSPLUS_MESH_BACKEND PROPERTY STRINGS MESHIO NATIVE KRATOS)
string(TOUPPER "${MESHIOPLUSPLUS_MESH_BACKEND}" _meshioplusplus_mesh_backend)

if(MESHIOPLUSPLUS_BUILD_PYTHON AND NOT _meshioplusplus_mesh_backend STREQUAL "MESHIO")
  message(FATAL_ERROR
    "meshio++: the pybind11 extension (MESHIOPLUSPLUS_BUILD_PYTHON=ON) requires "
    "MESHIOPLUSPLUS_MESH_BACKEND=MESHIO (got '${MESHIOPLUSPLUS_MESH_BACKEND}'). "
    "Configure with -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF to use the NATIVE/KRATOS backends.")
endif()

if(_meshioplusplus_mesh_backend STREQUAL "MESHIO")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_MESH_BACKEND_MESHIO)
elseif(_meshioplusplus_mesh_backend STREQUAL "NATIVE")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_MESH_BACKEND_NATIVE)
elseif(_meshioplusplus_mesh_backend STREQUAL "KRATOS")
  target_compile_definitions(meshioplusplus_core_obj PUBLIC MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)
else()
  message(FATAL_ERROR
    "meshio++: unknown MESHIOPLUSPLUS_MESH_BACKEND '${MESHIOPLUSPLUS_MESH_BACKEND}' "
    "(use MESHIO, NATIVE or KRATOS)")
endif()
message(STATUS "meshio++: mesh backend: ${_meshioplusplus_mesh_backend}")

if(MESHIOPLUSPLUS_BUILD_PYTHON)
  # The pybind11 extension: bindings + the shared core object library.
  file(GLOB MESHIOPLUSPLUS_BINDING_SOURCES bindings/python/*.cpp)
  pybind11_add_module(_core ${MESHIOPLUSPLUS_BINDING_SOURCES})
  target_include_directories(_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bindings/python)
  target_link_libraries(_core PRIVATE meshioplusplus_core_obj)

  # LDFLAGS env-var seeding of CMAKE_MODULE_LINKER_FLAGS was tried first and
  # silently didn't apply to this MODULE target, hence targeting it directly.
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
  file(GLOB MESHIOPLUSPLUS_JS_BINDING_SOURCES bindings/wasm/*.cpp)
  add_executable(meshioplusplus_wasm ${MESHIOPLUSPLUS_JS_BINDING_SOURCES})
  target_include_directories(meshioplusplus_wasm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bindings/wasm)
  target_link_libraries(meshioplusplus_wasm PRIVATE meshioplusplus_core_obj)
  # -sUSE_ZLIB=1 is not repeated here: meshioplusplus_core_obj already applies
  # it PUBLIC (both compile and link) when MESHIOPLUSPLUS_WITH_ZLIB is on, and
  # this target inherits it transitively via target_link_libraries above.
  # STACK_SIZE: Emscripten's default wasm stack is 64 KiB, which HDF5 and
  # netCDF-4 overrun -- their call chains are deep and their frames large (name
  # buffers, B-tree nodes, dataspace scratch). The overrun is SILENT: the stack
  # grows downward into the static data segment, so the first symptom is
  # unrelated global state going bad. Concretely, writing one Exodus file
  # clobbered libc++'s locale facets, and every subsequent `istream >> number`
  # in the module -- i.e. every ASCII reader, gmsh/obj/off/vtk included --
  # trapped with "null function or function signature mismatch". 4 MiB is
  # comfortably above what the deepest observed path needs (the minimum that
  # survived a netCDF-4 write was 128 KiB) and is only an address-space
  # reservation inside the growable heap.
  #
  # STACK_OVERFLOW_CHECK turns a future recurrence back into a loud abort
  # rather than silent corruption of whatever happens to sit below the stack.
  # It costs a canary check on stack allocation, which is the right trade for a
  # library whose whole job is not corrupting the user's data.
  target_link_options(
    meshioplusplus_wasm PRIVATE
    "--bind"
    "-sMODULARIZE=1"
    "-sEXPORT_ES6=1"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sSTACK_SIZE=4MB"
    "-sSTACK_OVERFLOW_CHECK=1"
    "-sFORCE_FILESYSTEM=1"
    "-sEXPORTED_RUNTIME_METHODS=['FS']")
  set_target_properties(meshioplusplus_wasm PROPERTIES OUTPUT_NAME "meshioplusplus_wasm"
                                                        SUFFIX ".mjs")
endif()

# --------------------------------------------------------------------------
# C API (bindings/c/): the installable `libmeshioplusplus` shared library +
# the pure-C header, the third flat binding over the same core (alongside
# WASM). Written against the uniform mesh API only, so it builds under every
# MESHIOPLUSPLUS_MESH_BACKEND. The object-library split lets the GoogleTest
# binary link the C API's objects directly (no RPATH/shared-lib coupling in
# the per-backend test legs); the shared lib is what gets installed/exported.
# --------------------------------------------------------------------------
option(MESHIOPLUSPLUS_BUILD_C_API "Build the installable libmeshioplusplus C API" OFF)
option(MESHIOPLUSPLUS_BUILD_FORTRAN "Build the Fortran module (implies MESHIOPLUSPLUS_BUILD_C_API)" OFF)
if(MESHIOPLUSPLUS_BUILD_FORTRAN AND NOT MESHIOPLUSPLUS_BUILD_C_API)
  message(STATUS "meshio++: MESHIOPLUSPLUS_BUILD_FORTRAN=ON force-enables MESHIOPLUSPLUS_BUILD_C_API")
  set(MESHIOPLUSPLUS_BUILD_C_API ON)
endif()

if(MESHIOPLUSPLUS_BUILD_C_API)
  include(GNUInstallDirs)  # before any use of CMAKE_INSTALL_*DIR below
  include(CMakePackageConfigHelpers)

  add_library(meshioplusplus_c_obj OBJECT bindings/c/c_api.cpp)
  set_target_properties(meshioplusplus_c_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_include_directories(meshioplusplus_c_obj
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bindings/c/include>)
  # MIO_SHARED/MIO_BUILDING drive the MIO_API dllexport/dllimport macro on
  # Windows; consumers of the installed library get MIO_SHARED from the
  # INTERFACE definition on the shared lib below.
  target_compile_definitions(meshioplusplus_c_obj
    PRIVATE MIO_BUILDING MIO_SHARED MIO_VERSION_STRING="${PROJECT_VERSION}")
  target_link_libraries(meshioplusplus_c_obj PUBLIC meshioplusplus_core_obj)

  # PRIVATE link: the object files are embedded, but the C++ core's usage
  # requirements (in-tree include dirs, backend macros, HDF5/... link deps)
  # stay out of the exported interface -- the installed surface is the C
  # header alone. Both object libraries must be linked DIRECTLY: CMake embeds
  # only directly-linked OBJECT libraries' objects (a transitive one, like
  # core_obj via c_obj, would contribute usage requirements but no objects,
  # leaving the .so with undefined registry/format symbols).
  add_library(meshioplusplus SHARED)
  target_link_libraries(meshioplusplus PRIVATE meshioplusplus_c_obj meshioplusplus_core_obj)
  target_include_directories(meshioplusplus
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/bindings/c/include>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  target_compile_definitions(meshioplusplus INTERFACE MIO_SHARED)
  set_target_properties(meshioplusplus PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 0)  # C ABI declared unstable pre-1.0 of the C API

  install(TARGETS meshioplusplus EXPORT meshioplusplusTargets
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES bindings/c/include/meshioplusplus/meshioplusplus.h
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/meshioplusplus)
  install(EXPORT meshioplusplusTargets NAMESPACE meshioplusplus::
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/meshioplusplus)
  configure_package_config_file(cmake/meshioplusplusConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/meshioplusplusConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/meshioplusplus)
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/meshioplusplusConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/meshioplusplusConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/meshioplusplusConfigVersion.cmake
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/meshioplusplus)

  # pkg-config, the lingua franca of HPC build systems. Libs.private carries
  # whatever optional deps this configure actually linked (static-link aid).
  set(MIO_PC_LIBS_PRIVATE "")
  if(HDF5_FOUND)
    string(APPEND MIO_PC_LIBS_PRIVATE " -lhdf5")
  endif()
  if(netCDF_FOUND OR NETCDF_LIBRARY)
    string(APPEND MIO_PC_LIBS_PRIVATE " -lnetcdf")
  endif()
  if(ZLIB_FOUND)
    string(APPEND MIO_PC_LIBS_PRIVATE " -lz")
  endif()
  if(MESHIOPLUSPLUS_WITH_KAHIP AND KAHIP_FOUND)
    string(APPEND MIO_PC_LIBS_PRIVATE " -lkahip")
  endif()
  configure_file(cmake/meshioplusplus.pc.in
                 ${CMAKE_CURRENT_BINARY_DIR}/meshioplusplus.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/meshioplusplus.pc
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
endif()

# --------------------------------------------------------------------------
# Fortran module (bindings/fortran/): a modern OO Fortran 2008 interface
# (`type(mio_mesh)` with type-bound procedures) layered on the C API via
# ISO_C_BINDING. The .f90 source is installed alongside the compiled .mod
# because .mod files are compiler-(major-version-)specific -- consumers on a
# different compiler recompile the module from source (the HDF5 approach).
# --------------------------------------------------------------------------
if(MESHIOPLUSPLUS_BUILD_FORTRAN)
  enable_language(Fortran)
  add_library(meshioplusplus_fortran SHARED bindings/fortran/meshioplusplus.f90)
  target_link_libraries(meshioplusplus_fortran PUBLIC meshioplusplus)
  set_target_properties(meshioplusplus_fortran PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 0
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules)
  target_include_directories(meshioplusplus_fortran INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran_modules>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/meshioplusplus/fortran>)

  install(TARGETS meshioplusplus_fortran EXPORT meshioplusplusTargets
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules/meshioplusplus.mod
                bindings/fortran/meshioplusplus.f90
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/meshioplusplus/fortran)
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
  # Keep `cmake --install` of a tests-enabled tree from installing gtest
  # next to the C API (relevant since MESHIOPLUSPLUS_BUILD_C_API).
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  file(GLOB MESHIOPLUSPLUS_TEST_SOURCES tests/cpp/*.cpp)
  if(NOT MESHIOPLUSPLUS_BUILD_C_API)
    list(REMOVE_ITEM MESHIOPLUSPLUS_TEST_SOURCES
         ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/test_c_api.cpp)
  endif()
  add_executable(meshioplusplus_tests ${MESHIOPLUSPLUS_TEST_SOURCES})
  target_include_directories(meshioplusplus_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp)
  # The Polyscope payload builder is CLI-side but deliberately Polyscope-free,
  # so its invariants are tested here in the default build -- no OpenGL needed.
  target_sources(meshioplusplus_tests PRIVATE src/cpp/cli/view_payload.cpp)
  target_link_libraries(meshioplusplus_tests PRIVATE meshioplusplus_core_obj GTest::gtest_main)
  if(MESHIOPLUSPLUS_BUILD_C_API)
    # Link the C API's objects directly (not the shared lib): the gtest binary
    # then needs no install RPATH and the per-backend CI legs stay one-target.
    target_link_libraries(meshioplusplus_tests PRIVATE meshioplusplus_c_obj)
  endif()

  include(GoogleTest)
  gtest_discover_tests(meshioplusplus_tests)

  if(MESHIOPLUSPLUS_BUILD_FORTRAN)
    # Plain Fortran program, nonzero exit on any failed check; argv[1] is a
    # path prefix for the files it writes.
    add_executable(meshioplusplus_fortran_test tests/fortran/test_fortran_api.f90)
    target_link_libraries(meshioplusplus_fortran_test PRIVATE meshioplusplus_fortran)
    add_test(NAME fortran_api
             COMMAND meshioplusplus_fortran_test ${CMAKE_CURRENT_BINARY_DIR}/fortran_test_out)
  endif()
endif()

# --------------------------------------------------------------------------
# C++ backend benchmark (src/cpp/benchmark/bench_backends.cpp): one binary per
# mesh backend (the backend is a compile-time choice); benchmark/
# bench_backends.sh builds all three variants and collates the CSV output.
# No external benchmark framework -- std::chrono, warmup + median-of-N.
# --------------------------------------------------------------------------
option(MESHIOPLUSPLUS_BUILD_BENCHMARKS "Build the C++ mesh-backend benchmark binary" OFF)
if(MESHIOPLUSPLUS_BUILD_BENCHMARKS)
  add_executable(meshioplusplus_bench src/cpp/benchmark/bench_backends.cpp)
  target_link_libraries(meshioplusplus_bench PRIVATE meshioplusplus_core_obj)
endif()

# --------------------------------------------------------------------------
# Native command-line binary (src/cpp/cli/main.cpp): a Python-free CLI over
# meshioplusplus_core_obj, mirroring the Python console script's verbs
# (convert/info/ascii/binary/compress/decompress/quality/extract-surface/
# reorder/diff). Installed executable name is `meshioplusplus`. Opt-in like the
# other executables; works under any mesh/parallel backend and optional-dep set.
# --------------------------------------------------------------------------
option(MESHIOPLUSPLUS_BUILD_CLI "Build the native command-line binary" OFF)
if(MESHIOPLUSPLUS_BUILD_CLI)
  add_executable(meshioplusplus_cli src/cpp/cli/main.cpp)
  target_sources(meshioplusplus_cli PRIVATE src/cpp/cli/view_payload.cpp)
  target_link_libraries(meshioplusplus_cli PRIVATE meshioplusplus_core_obj)
  target_compile_definitions(meshioplusplus_cli
    PRIVATE MESHIOPLUSPLUS_CLI_VERSION="${PROJECT_VERSION}")

  # The viewer half. PRIVATE on purpose -- nothing else in the project may
  # acquire an OpenGL dependency through it. EXCLUDE_FROM_ALL keeps Polyscope's
  # own examples and tests out of the build.
  if(MESHIOPLUSPLUS_WITH_POLYSCOPE AND NOT EMSCRIPTEN)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/third_party/polyscope/CMakeLists.txt")
      message(STATUS "meshio++: Polyscope viewer enabled (CLI only)")
      add_subdirectory(src/cpp/third_party/polyscope EXCLUDE_FROM_ALL)
      target_sources(meshioplusplus_cli PRIVATE src/cpp/cli/polyscope_view.cpp)
      target_link_libraries(meshioplusplus_cli PRIVATE polyscope)
      target_compile_definitions(meshioplusplus_cli PRIVATE MESHIOPLUSPLUS_HAS_POLYSCOPE)
    else()
      message(STATUS
        "meshio++: Polyscope requested but src/cpp/third_party/polyscope is empty - "
        "run `git submodule update --init --recursive`; view/screenshot will raise")
    endif()
  endif()
  set_target_properties(meshioplusplus_cli PROPERTIES OUTPUT_NAME meshioplusplus)
  if(MESHIOPLUSPLUS_STATIC_RUNTIME)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_link_options(meshioplusplus_cli PRIVATE -static-libgcc -static-libstdc++)
    elseif(MSVC)
      # CMP0091 (NEW via the top-of-file cmake_minimum_required range) makes this
      # property authoritative regardless of when the target was created: static
      # CRT (/MT, /MTd for Debug) so the binary needs no vcredist at runtime.
      set_property(TARGET meshioplusplus_cli PROPERTY
        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
  endif()
  include(GNUInstallDirs)
  install(TARGETS meshioplusplus_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
