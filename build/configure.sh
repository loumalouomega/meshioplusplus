#!/bin/sh
# configure.sh — configure (and optionally build) the meshio++ C++ core
# standalone on Linux/macOS. Run from anywhere; build trees are created next
# to this script (build/cpp-<build-type>).
#
#   ./configure.sh                          # defaults: STL backend, Release,
#                                           # HDF5/netCDF/zlib auto-detected
#   ./configure.sh --backend OPENMP --tests --build
#   ./configure.sh --backend TBB --tbb-dir /opt/intel/oneapi/tbb/latest/lib/cmake/tbb
#   ./configure.sh --mesh-backend NATIVE --tests --build   # standalone, no Python
#
# The equivalent Python-package install is printed at the end.

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_DIR=$(dirname -- "$SCRIPT_DIR")

BACKEND="STL"
MESH_BACKEND="MESHIO"
BUILD_TYPE="Release"
WITH_HDF5="ON"
WITH_NETCDF="ON"
WITH_ZLIB="ON"
WITH_GIDPOST="ON"
WITH_POLYSCOPE="OFF"
TESTS="OFF"
C_API="OFF"
FORTRAN="OFF"
INSTALL_CPP="OFF"
CPP_BACKENDS=""
CLI="OFF"
DO_BUILD="no"
PYTHON_EXE=""
TBB_DIR=""

usage() {
    cat <<EOF
Usage: $0 [options]
  --backend <SEQ|STL|OPENMP|TBB|KOKKOS>  parallel backend (default: STL)
  --mesh-backend <MESHIO|NATIVE|KRATOS>
                                  in-memory mesh backend (default: MESHIO).
                                  NATIVE/KRATOS imply -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF
                                  (the pybind11 extension requires MESHIO)
  --build-type <type>             CMake build type (default: Release)
  --with-hdf5 / --without-hdf5    HDF5-backed formats (default: on, auto-detected)
  --with-netcdf / --without-netcdf
  --with-zlib / --without-zlib
  --with-gidpost / --without-gidpost
                                  GiD postprocess writer (default: on; needs zlib)
  --with-polyscope                native viewer for the CLI's view/screenshot
                                  (default: off; needs OpenGL/GLFW and
                                  git submodule update --init --recursive)
  --tests                         also build the GoogleTest suite (CTest)
  --c-api                         build the installable libmeshioplusplus C API
  --fortran                       build the Fortran module (implies --c-api)
  --install-cpp                   build the installable C++ API
                                  (meshioplusplus::core*, see doc/cpp_api.md)
  --cpp-backends <LIST>           mesh backends to build as
                                  meshioplusplus::core_<backend> (comma- or
                                  semicolon-separated; default MESHIO,NATIVE,KRATOS;
                                  implies --install-cpp)
  --cli                           build the native command-line binary (meshioplusplus)
  --build                         run the build after configuring
  --python <exe>                  Python executable (default: auto)
  --tbb-dir <path>                TBBConfig.cmake dir (e.g. oneAPI)
  -h, --help                      this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --backend) BACKEND="$2"; shift 2 ;;
        --mesh-backend) MESH_BACKEND=$(echo "$2" | tr '[:lower:]' '[:upper:]'); shift 2 ;;
        --build-type) BUILD_TYPE="$2"; shift 2 ;;
        --with-hdf5) WITH_HDF5="ON"; shift ;;
        --without-hdf5) WITH_HDF5="OFF"; shift ;;
        --with-netcdf) WITH_NETCDF="ON"; shift ;;
        --without-netcdf) WITH_NETCDF="OFF"; shift ;;
        --with-zlib) WITH_ZLIB="ON"; shift ;;
        --with-polyscope) WITH_POLYSCOPE="ON"; shift ;;
        --without-polyscope) WITH_POLYSCOPE="OFF"; shift ;;
        --without-zlib) WITH_ZLIB="OFF"; shift ;;
        --with-gidpost) WITH_GIDPOST="ON"; shift ;;
        --without-gidpost) WITH_GIDPOST="OFF"; shift ;;
        --tests) TESTS="ON"; shift ;;
        --c-api) C_API="ON"; shift ;;
        --fortran) FORTRAN="ON"; C_API="ON"; shift ;;
        --install-cpp) INSTALL_CPP="ON"; shift ;;
        --cpp-backends) CPP_BACKENDS=$(echo "$2" | tr '[:lower:],' '[:upper:];'); INSTALL_CPP="ON"; shift 2 ;;
        --cli) CLI="ON"; shift ;;
        --build) DO_BUILD="yes"; shift ;;
        --python) PYTHON_EXE="$2"; shift 2 ;;
        --tbb-dir) TBB_DIR="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# gidpostInt.h includes <zlib.h> unconditionally and the binary flavour is
# always deflated, so the GiD writer cannot be built without zlib.
if [ "$WITH_GIDPOST" = "ON" ] && [ "$WITH_ZLIB" = "OFF" ]; then
    echo "error: --with-gidpost needs zlib (gidpost deflates unconditionally)." >&2
    echo "       Use --without-gidpost as well, or drop --without-zlib." >&2
    exit 1
fi

# Python: prefer an in-repo venv, then python3.
if [ -z "$PYTHON_EXE" ]; then
    if [ -x "$SOURCE_DIR/.venv/bin/python" ]; then
        PYTHON_EXE="$SOURCE_DIR/.venv/bin/python"
    else
        PYTHON_EXE=$(command -v python3 || command -v python)
    fi
fi

# Non-MESHIO backends get their own build tree (and no Python extension).
BUILD_PYTHON="ON"
TREE_SUFFIX=""
if [ "$MESH_BACKEND" != "MESHIO" ]; then
    BUILD_PYTHON="OFF"
    TREE_SUFFIX="-$(echo "$MESH_BACKEND" | tr '[:upper:]' '[:lower:]')"
fi

BUILD_DIR="$SCRIPT_DIR/cpp-$(echo "$BUILD_TYPE" | tr '[:upper:]' '[:lower:]')$TREE_SUFFIX"

GENERATOR=""
if command -v ninja >/dev/null 2>&1; then
    GENERATOR="-G Ninja"
fi

set -- \
    -S "$SOURCE_DIR" -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DMESHIOPLUSPLUS_PARALLEL_BACKEND="$BACKEND" \
    -DMESHIOPLUSPLUS_MESH_BACKEND="$MESH_BACKEND" \
    -DMESHIOPLUSPLUS_BUILD_PYTHON="$BUILD_PYTHON" \
    -DMESHIOPLUSPLUS_WITH_HDF5="$WITH_HDF5" \
    -DMESHIOPLUSPLUS_WITH_NETCDF="$WITH_NETCDF" \
    -DMESHIOPLUSPLUS_WITH_ZLIB="$WITH_ZLIB" \
    -DMESHIOPLUSPLUS_WITH_GIDPOST="$WITH_GIDPOST" \
        -DMESHIOPLUSPLUS_WITH_POLYSCOPE="$WITH_POLYSCOPE" \
    -DMESHIOPLUSPLUS_BUILD_TESTS="$TESTS" \
    -DMESHIOPLUSPLUS_BUILD_C_API="$C_API" \
    -DMESHIOPLUSPLUS_INSTALL_CPP="$INSTALL_CPP" \
    -DMESHIOPLUSPLUS_BUILD_CLI="$CLI" \
    -DMESHIOPLUSPLUS_BUILD_FORTRAN="$FORTRAN" \
    -DPython_EXECUTABLE="$PYTHON_EXE"

# Only pass the backend list when the user gave one, so the CMake cache
# default (MESHIO;NATIVE;KRATOS) stays the single source of truth.
if [ -n "$CPP_BACKENDS" ]; then
    set -- "$@" -DMESHIOPLUSPLUS_INSTALL_CPP_BACKENDS="$CPP_BACKENDS"
fi

# pybind11 from the chosen Python, when available.
PYBIND11_DIR=$("$PYTHON_EXE" -c "import pybind11; print(pybind11.get_cmake_dir())" 2>/dev/null || true)
if [ -n "$PYBIND11_DIR" ]; then
    set -- "$@" -Dpybind11_DIR="$PYBIND11_DIR"
fi
if [ -n "$TBB_DIR" ]; then
    set -- "$@" -DTBB_DIR="$TBB_DIR"
fi

echo "== meshio++ configure =="
echo "  source:    $SOURCE_DIR"
echo "  build:     $BUILD_DIR"
echo "  type:      $BUILD_TYPE"
echo "  backend:   $BACKEND"
echo "  mesh:      $MESH_BACKEND (Python extension: $BUILD_PYTHON)"
echo "  HDF5:      $WITH_HDF5   netCDF: $WITH_NETCDF   zlib: $WITH_ZLIB"
echo "  gidpost:   $WITH_GIDPOST  (GiD postprocess writer)"
echo "  Polyscope: $WITH_POLYSCOPE  (CLI viewer)"
echo "  tests:     $TESTS"
echo "  C API:     $C_API   Fortran: $FORTRAN"
echo "  C++ API:   $INSTALL_CPP${CPP_BACKENDS:+  (backends: $CPP_BACKENDS)}"
echo "  CLI:       $CLI"
echo "  python:    $PYTHON_EXE"
echo

# shellcheck disable=SC2086
cmake $GENERATOR "$@"

echo
echo "== next steps =="
echo "  cmake --build \"$BUILD_DIR\" -j"
if [ "$TESTS" = "ON" ]; then
    echo "  ctest --test-dir \"$BUILD_DIR\" --output-on-failure"
fi
if [ "$C_API" = "ON" ] && [ "$INSTALL_CPP" = "ON" ]; then
    echo "  cmake --install \"$BUILD_DIR\" --prefix <prefix>   # C + C++ APIs (see doc/c_api.md, doc/cpp_api.md)"
elif [ "$C_API" = "ON" ]; then
    echo "  cmake --install \"$BUILD_DIR\" --prefix <prefix>   # libmeshioplusplus + headers (see doc/c_api.md)"
elif [ "$INSTALL_CPP" = "ON" ]; then
    echo "  cmake --install \"$BUILD_DIR\" --prefix <prefix>   # C++ API: include/ lib/ cmake config (see doc/cpp_api.md)"
fi
if [ "$CLI" = "ON" ]; then
    echo "  \"$BUILD_DIR/meshioplusplus\" --help   # native CLI (no Python)"
fi
echo
if [ "$MESH_BACKEND" = "MESHIO" ]; then
    echo "Python package (editable) with the same options:"
    echo "  CMAKE_ARGS=\"-DMESHIOPLUSPLUS_PARALLEL_BACKEND=$BACKEND -DMESHIOPLUSPLUS_WITH_HDF5=$WITH_HDF5 -DMESHIOPLUSPLUS_WITH_NETCDF=$WITH_NETCDF -DMESHIOPLUSPLUS_WITH_ZLIB=$WITH_ZLIB\" \\"
    echo "    pip install --no-build-isolation -e \"$SOURCE_DIR\""
else
    echo "Note: the $MESH_BACKEND mesh backend is standalone-C++ only; the Python"
    echo "package always uses MESHIO (PyPI wheels are unaffected)."
fi

if [ "$DO_BUILD" = "yes" ]; then
    echo
    echo "== building =="
    cmake --build "$BUILD_DIR" -j
fi
