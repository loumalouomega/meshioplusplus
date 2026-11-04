#!/bin/sh
# build-cutcells.sh -- build the pinned CutCells (MIT) for the optional
# level-set cutting backend (-DMESHIOPLUSPLUS_WITH_CUTCELLS=ON, doc/cut.md).
#
#   build/build-cutcells.sh PREFIX            # static, position-independent
#   build/build-cutcells.sh PREFIX --shared   # a shared libcutcells
#
# CutCells is cloned at CUTCELLS_COMMIT (build/deps-pins.sh) and built without
# OpenMP, so a static meshio++ (the wheels, the CLI binaries) carries no OpenMP
# runtime of its own; meshio++ calls it serially.
#
# Then configure meshio++ with -DMESHIOPLUSPLUS_WITH_CUTCELLS=ON
# -DCUTCELLS_ROOT=PREFIX (or build/configure.sh --with-cutcells --cutcells-root PREFIX).
#
# Like build-mmg.sh, a shell script and not CMake: meshio++'s CMake never
# downloads anything, it only finds an already-built dependency.

set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$HERE/deps-pins.sh"

PREFIX=""
SHARED="OFF"
FORCE="no"
JOBS=""
WORK=""

usage() {
    cat <<USAGE
Usage: build/build-cutcells.sh PREFIX [options]
  --shared         build a shared library (default: static, position-independent)
  --work DIR       clone and build directory (default: PREFIX/_work)
  --jobs N         parallel build jobs (default: all cores)
  --force          rebuild even when PREFIX already holds CutCells
Builds CutCells $CUTCELLS_VERSION at $CUTCELLS_COMMIT, without OpenMP.
USAGE
}

while [ $# -gt 0 ]; do
    case "$1" in
        --shared) SHARED="ON"; shift ;;
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
[ -n "$JOBS" ] || JOBS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo "${NUMBER_OF_PROCESSORS:-2}")

if [ "$FORCE" = "no" ] && [ -f "$PREFIX/lib/cmake/cutcells/CutCellsConfig.cmake" ]; then
    echo "build-cutcells: CutCells is already in $PREFIX (--force rebuilds)"
    exit 0
fi

PLATFORM_CMAKE=""
case "$(uname -s)" in
    Darwin) PLATFORM_CMAKE="-DCMAKE_OSX_ARCHITECTURES=x86_64;arm64 -DCMAKE_OSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-13.3}" ;;
esac

SRC="$WORK/cutcells"
if [ ! -d "$SRC/.git" ]; then
    git clone --quiet https://github.com/sclaus2/CutCells.git "$SRC"
fi
git -C "$SRC" fetch --quiet origin "$CUTCELLS_COMMIT" 2>/dev/null || true
git -C "$SRC" checkout --quiet "$CUTCELLS_COMMIT"
[ "$(git -C "$SRC" rev-parse HEAD)" = "$CUTCELLS_COMMIT" ] || {
    echo "build-cutcells: the clone is not at $CUTCELLS_COMMIT" >&2; exit 1; }

# MSVC has no __int128, which CutCells' quadrays code uses for exact rational
# sums: the patch gives it a small 128-bit integer there and leaves every other
# compiler on the original typedef. Applied once (a rerun finds it applied).
PATCH="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)/patches/cutcells-$CUTCELLS_VERSION-msvc-int128.patch"
if git -C "$SRC" apply --reverse --check "$PATCH" 2>/dev/null; then
    :   # already applied
else
    git -C "$SRC" apply "$PATCH"
fi

# shellcheck disable=SC2086  # PLATFORM_CMAKE is meant to word-split
cmake -S "$SRC/cpp" -B "$WORK/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_INSTALL_LIBDIR=lib \
    -DBUILD_SHARED_LIBS="$SHARED" -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON -DBUILD_TESTING=OFF $PLATFORM_CMAKE
cmake --build "$WORK/build" --config Release --parallel "$JOBS"
cmake --install "$WORK/build" --config Release
echo "build-cutcells: CutCells $CUTCELLS_VERSION installed in $PREFIX"
