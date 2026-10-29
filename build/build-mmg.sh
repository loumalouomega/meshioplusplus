#!/bin/sh
# build-mmg.sh -- build the pinned MMG (mmg2d, mmgs, mmg3d) as shared libraries
# for the optional MMG backend (-DMESHIOPLUSPLUS_WITH_MMG=ON, roadmap 3.2).
# MMG is LGPL-3.0: meshio++ links it dynamically only, never vendors it, and
# leaves it out of the release wheels and CLI binaries (doc/mmg.md).
#
#   build/build-mmg.sh PREFIX                # MMG alone: adapt, discretize_levelset, mesh_domain
#   build/build-mmg.sh PREFIX --with-elas    # also ELAS, for move (GPL-3.0, see below)
#
# libmmg is given an RPATH to PREFIX/lib, so it finds ELAS there by itself.
#
# Then configure meshio++ with -DMESHIOPLUSPLUS_WITH_MMG=ON -Dmmg_DIR=PREFIX/lib/cmake/mmg
# (or build/configure.sh --with-mmg --mmg-root PREFIX), and put PREFIX/lib on
# LD_LIBRARY_PATH when running a Python build, whose RPATH scikit-build-core strips.
#
# --with-elas builds the ISCD LinearElasticity library and its Commons first
# and builds MMG against them, which is what MMG3D_mmg3dmov and MMG2D_mmg2dmov
# need. Both are GPL-3.0, so an MMG built with them -- and a meshio++ linked to
# that MMG -- is a GPL-3.0 combination when it is redistributed.
#
# Like build-occt.sh, a shell script and not CMake: meshio++'s CMake never
# downloads anything, it only finds an already-built dependency.

set -eu

MMG_VERSION="v5.8.0"
COMMONS_COMMIT="61c211a4717fb36621dd36700ffc1c9b8c5bb93c"
ELAS_COMMIT="58ace131fd1af9293ded9f84f2a6447147a8ecdc"

PREFIX=""
WITH_ELAS="no"
FORCE="no"
JOBS=""
WORK=""

usage() {
    cat <<EOF
Usage: build/build-mmg.sh PREFIX [options]
  --with-elas      also build ELAS (and its Commons) and MMG against it, for move (GPL-3.0)
  --work DIR       download and build directory (default: PREFIX/_work)
  --jobs N         parallel build jobs (default: all cores)
  --force          rebuild even when PREFIX already holds the libraries
Builds MMG $MMG_VERSION shared, 32-bit MMG5_int, without SCOTCH or VTK.
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --with-elas) WITH_ELAS="yes"; shift ;;
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
[ -n "$JOBS" ] || JOBS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)

if [ "$FORCE" = "no" ] && [ -f "$PREFIX/lib/cmake/mmg/mmgConfig.cmake" ]; then
    if [ "$WITH_ELAS" = "no" ] || [ -n "$(ls "$PREFIX"/lib/libElas.* 2>/dev/null)" ]; then
        echo "build-mmg: $PREFIX already holds MMG (--force to rebuild)"
        exit 0
    fi
fi

fetch() {  # fetch URL REF DIR: a shallow checkout of one tag or commit
    if [ ! -d "$3/.git" ]; then
        git init -q "$3"
        git -C "$3" remote add origin "$1"
    fi
    git -C "$3" fetch -q --depth 1 origin "$2"
    git -C "$3" checkout -q FETCH_HEAD
}

MMG_CMAKE_ELAS="-DUSE_ELAS=OFF"
if [ "$WITH_ELAS" = "yes" ]; then
    # ISCD's CMake installs to $HOME/{lib,include}: point HOME at the prefix.
    fetch https://github.com/ISCDtoolbox/Commons "$COMMONS_COMMIT" "$WORK/commons"
    fetch https://github.com/ISCDtoolbox/LinearElasticity "$ELAS_COMMIT" "$WORK/elas"
    for dep in commons elas; do
        HOME="$PREFIX" cmake -S "$WORK/$dep" -B "$WORK/$dep-build" -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_POLICY_VERSION_MINIMUM=3.5
        HOME="$PREFIX" cmake --build "$WORK/$dep-build" -j "$JOBS"
        HOME="$PREFIX" cmake --install "$WORK/$dep-build"
    done
    MMG_CMAKE_ELAS="-DUSE_ELAS=ON -DELAS_DIR=$PREFIX"
fi

fetch https://github.com/MmgTools/mmg "$MMG_VERSION" "$WORK/mmg"
# shellcheck disable=SC2086
cmake -S "$WORK/mmg" -B "$WORK/mmg-build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_INSTALL_RPATH="$PREFIX/lib" \
    -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF \
    -DUSE_SCOTCH=OFF -DUSE_VTK=OFF -DMMG5_INT=int32_t $MMG_CMAKE_ELAS
cmake --build "$WORK/mmg-build" -j "$JOBS"
cmake --install "$WORK/mmg-build"
echo "build-mmg: MMG $MMG_VERSION installed in $PREFIX (ELAS: $WITH_ELAS)"
