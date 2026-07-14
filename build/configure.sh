#!/bin/sh
# configure.sh — configure (and optionally build) the meshio C++ core
# standalone on Linux/macOS. Run from anywhere; build trees are created next
# to this script (build/cpp-<build-type>).
#
#   ./configure.sh                          # defaults: STL backend, Release,
#                                           # HDF5/netCDF/zlib auto-detected
#   ./configure.sh --backend OPENMP --tests --build
#   ./configure.sh --backend TBB --tbb-dir /opt/intel/oneapi/tbb/latest/lib/cmake/tbb
#
# The equivalent Python-package install is printed at the end.

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_DIR=$(dirname -- "$SCRIPT_DIR")

BACKEND="STL"
BUILD_TYPE="Release"
WITH_HDF5="ON"
WITH_NETCDF="ON"
WITH_ZLIB="ON"
TESTS="OFF"
DO_BUILD="no"
PYTHON_EXE=""
TBB_DIR=""

usage() {
    cat <<EOF
Usage: $0 [options]
  --backend <SEQ|STL|OPENMP|TBB>  parallel backend (default: STL)
  --build-type <type>             CMake build type (default: Release)
  --with-hdf5 / --without-hdf5    HDF5-backed formats (default: on, auto-detected)
  --with-netcdf / --without-netcdf
  --with-zlib / --without-zlib
  --tests                         also build the GoogleTest suite (CTest)
  --build                         run the build after configuring
  --python <exe>                  Python executable (default: auto)
  --tbb-dir <path>                TBBConfig.cmake dir (e.g. oneAPI)
  -h, --help                      this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --backend) BACKEND="$2"; shift 2 ;;
        --build-type) BUILD_TYPE="$2"; shift 2 ;;
        --with-hdf5) WITH_HDF5="ON"; shift ;;
        --without-hdf5) WITH_HDF5="OFF"; shift ;;
        --with-netcdf) WITH_NETCDF="ON"; shift ;;
        --without-netcdf) WITH_NETCDF="OFF"; shift ;;
        --with-zlib) WITH_ZLIB="ON"; shift ;;
        --without-zlib) WITH_ZLIB="OFF"; shift ;;
        --tests) TESTS="ON"; shift ;;
        --build) DO_BUILD="yes"; shift ;;
        --python) PYTHON_EXE="$2"; shift 2 ;;
        --tbb-dir) TBB_DIR="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# Python: prefer an in-repo venv, then python3.
if [ -z "$PYTHON_EXE" ]; then
    if [ -x "$SOURCE_DIR/.venv/bin/python" ]; then
        PYTHON_EXE="$SOURCE_DIR/.venv/bin/python"
    else
        PYTHON_EXE=$(command -v python3 || command -v python)
    fi
fi

BUILD_DIR="$SCRIPT_DIR/cpp-$(echo "$BUILD_TYPE" | tr '[:upper:]' '[:lower:]')"

GENERATOR=""
if command -v ninja >/dev/null 2>&1; then
    GENERATOR="-G Ninja"
fi

set -- \
    -S "$SOURCE_DIR" -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DMESHIO_PARALLEL_BACKEND="$BACKEND" \
    -DMESHIO_WITH_HDF5="$WITH_HDF5" \
    -DMESHIO_WITH_NETCDF="$WITH_NETCDF" \
    -DMESHIO_WITH_ZLIB="$WITH_ZLIB" \
    -DMESHIO_BUILD_TESTS="$TESTS" \
    -DPython_EXECUTABLE="$PYTHON_EXE"

# pybind11 from the chosen Python, when available.
PYBIND11_DIR=$("$PYTHON_EXE" -c "import pybind11; print(pybind11.get_cmake_dir())" 2>/dev/null || true)
if [ -n "$PYBIND11_DIR" ]; then
    set -- "$@" -Dpybind11_DIR="$PYBIND11_DIR"
fi
if [ -n "$TBB_DIR" ]; then
    set -- "$@" -DTBB_DIR="$TBB_DIR"
fi

echo "== meshio configure =="
echo "  source:    $SOURCE_DIR"
echo "  build:     $BUILD_DIR"
echo "  type:      $BUILD_TYPE"
echo "  backend:   $BACKEND"
echo "  HDF5:      $WITH_HDF5   netCDF: $WITH_NETCDF   zlib: $WITH_ZLIB"
echo "  tests:     $TESTS"
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
echo
echo "Python package (editable) with the same options:"
echo "  CMAKE_ARGS=\"-DMESHIO_PARALLEL_BACKEND=$BACKEND -DMESHIO_WITH_HDF5=$WITH_HDF5 -DMESHIO_WITH_NETCDF=$WITH_NETCDF -DMESHIO_WITH_ZLIB=$WITH_ZLIB\" \\"
echo "    pip install --no-build-isolation -e \"$SOURCE_DIR\""

if [ "$DO_BUILD" = "yes" ]; then
    echo
    echo "== building =="
    cmake --build "$BUILD_DIR" -j
fi
