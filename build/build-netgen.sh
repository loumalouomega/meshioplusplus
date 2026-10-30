#!/bin/sh
# build-netgen.sh -- build the pinned Netgen (nglib and ngcore) as shared
# libraries for the Netgen tetrahedralize backend (-DMESHIOPLUSPLUS_WITH_NETGEN=ON
# in a source build, the meshioplusplus-netgen add-on wheel, roadmap 4.7.1).
# Netgen is LGPL-2.1: meshio++ links it dynamically only and never vendors it.
#
#   build/build-netgen.sh PREFIX
#
# Then configure meshio++ with -DMESHIOPLUSPLUS_WITH_NETGEN=ON -DNETGEN_ROOT=PREFIX
# -DMESHIOPLUSPLUS_NETGEN_VERSION=6.2.2401 (cmake/FindNetgen.cmake reads
# NETGEN_ROOT; nglib.h carries no version, so the build is told it), or build the
# add-on with NETGEN_ROOT set and that definition. The pin is the version the
# native CI job and the full WebAssembly artifact use (build-wasm-deps.sh).
#
# Like build-occt.sh, a shell script and not CMake: meshio++'s CMake never
# downloads anything, it only finds an already-built dependency.

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$SCRIPT_DIR/deps-pins.sh"   # NETGEN_VERSION, NETGEN_COMMIT

PREFIX=""
FORCE="no"
JOBS=""
WORK=""

usage() {
    cat <<EOF
Usage: build/build-netgen.sh PREFIX [options]
  --work DIR       download and build directory (default: PREFIX/_work)
  --jobs N         parallel build jobs (default: all cores)
  --force          rebuild even when PREFIX already holds the libraries
Builds Netgen $NETGEN_VERSION nglib and ngcore shared, without the GUI, Python,
OpenCASCADE or MPI.
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --work) WORK="$2"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        --force) FORCE="yes"; shift ;;
        -h|--help) usage; exit 0 ;;
        -*) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
        *) [ -z "$PREFIX" ] || { echo "one PREFIX only" >&2; exit 2; }; PREFIX="$1"; shift ;;
    esac
done
[ -n "$PREFIX" ] || { usage >&2; exit 2; }
mkdir -p "$PREFIX"
PREFIX=$(CDPATH= cd -- "$PREFIX" && pwd)
[ -n "$WORK" ] || WORK="$PREFIX/_work"
mkdir -p "$WORK"
[ -n "$JOBS" ] || JOBS=$( (command -v nproc >/dev/null 2>&1 && nproc) ||
                          (sysctl -n hw.ncpu 2>/dev/null) || echo 2 )

if [ "$FORCE" = "no" ] && { [ -f "$PREFIX/lib/libnglib.so" ] || [ -f "$PREFIX/lib/libnglib.dylib" ] ||
                           [ -f "$PREFIX/lib/nglib.lib" ] || [ -f "$PREFIX/bin/nglib.dll" ]; }; then
    echo "== Netgen $NETGEN_VERSION already installed in $PREFIX, skipping =="
    exit 0
fi

echo "== Netgen $NETGEN_VERSION (shared) into $PREFIX =="
SRC="$WORK/netgen-$NETGEN_VERSION"
if [ ! -d "$SRC/.git" ]; then
    rm -rf "$SRC"
    git -c advice.detachedHead=false clone --quiet --depth 1 --branch "v$NETGEN_VERSION" \
        https://github.com/NGSolve/netgen.git "$SRC"
fi
got=$(git -C "$SRC" rev-parse HEAD)
if [ "$got" != "$NETGEN_COMMIT" ]; then
    echo "error: netgen v$NETGEN_VERSION is $got, expected $NETGEN_COMMIT" >&2
    exit 1
fi

# nglib needs none of Netgen's optional components, and its version comes in by
# option because the build is not from a describe-able checkout.
cmake -S "$SRC" -B "$WORK/build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DNETGEN_VERSION_GIT="v$NETGEN_VERSION-0-g$(echo "$NETGEN_COMMIT" | cut -c1-7)" \
    -DUSE_SUPERBUILD=OFF -DUSE_GUI=OFF -DUSE_PYTHON=OFF -DUSE_OCC=OFF \
    -DUSE_MPI=OFF -DUSE_NATIVE_ARCH=OFF -DUSE_JPEG=OFF -DUSE_MPEG=OFF \
    -DUSE_CGNS=OFF -DBUILD_STUB_FILES=OFF -DENABLE_UNIT_TESTS=OFF \
    -DNGLIB_LIBRARY_TYPE=SHARED -DNGCORE_LIBRARY_TYPE=SHARED
cmake --build "$WORK/build" --config Release -j "$JOBS"
cmake --install "$WORK/build" --config Release
rm -rf "$WORK/build" "$PREFIX/share/netgen"
echo "== done: -DNETGEN_ROOT=$PREFIX -DMESHIOPLUSPLUS_NETGEN_VERSION=$NETGEN_VERSION =="
