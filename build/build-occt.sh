#!/bin/sh
# build-occt.sh -- build the pinned OpenCASCADE (build/occt-recipe.sh) as a
# static, position-independent install the core links into one binary. This is
# how the PyPI wheels and the release CLI binaries carry the kernel: linked in,
# so no wheel repair tool, RPATH or second C++ runtime is involved, and the
# LGPL obligations are met by the notices and the relink instructions in
# THIRD_PARTY_NOTICES.md (this script is the relink recipe).
#
#   build/build-occt.sh PREFIX                      # Linux, macOS, Windows (Git Bash)
#   build/build-occt.sh PREFIX --universal          # macOS x86_64 + arm64, 13.3
#   build/build-occt.sh PREFIX --msvc-static-runtime  # /MT, for the CLI binaries
#   build/build-occt.sh PREFIX --shared             # shared libraries, for the
#                                                   # meshioplusplus-occt add-on wheel
#
# The add-on ships OpenCASCADE as shared objects beside its own extension, so a
# recipient can swap the LGPL library without rebuilding anything.
#
# Then configure meshio++ with -DOpenCASCADE_DIR=PREFIX/lib/cmake/opencascade
# (MESHIOPLUSPLUS_WITH_OCCT is on by default). A source build can equally use a
# system OpenCASCADE 7.6 or later, linked dynamically, instead.
#
# Like build-wasm-deps.sh, a shell script and not CMake: meshio++'s CMake never
# downloads anything, it only finds an already-built dependency.

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OCCT_RECIPE_DIR="$SCRIPT_DIR"
. "$SCRIPT_DIR/occt-recipe.sh"

PREFIX=""
UNIVERSAL="no"
STATIC_CRT="no"
SHARED="no"
FORCE="no"
JOBS=""
WORK=""

usage() {
    cat <<EOF
Usage: build/build-occt.sh PREFIX [options]
  --universal             macOS: build x86_64 and arm64 (deployment target 13.3)
  --msvc-static-runtime   MSVC: link the static C runtime (/MT), as the CLI binaries do
  --shared                build shared libraries (the meshioplusplus-occt add-on wheel)
  --work DIR              download and build directory (default: PREFIX/_work)
  --jobs N                parallel build jobs (default: all cores)
  --force                 rebuild even when PREFIX already holds the libraries
Builds OpenCASCADE $OCCT_VERSION ($OCCT_TOOLKITS and their closure), static unless --shared.
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --universal) UNIVERSAL="yes"; shift ;;
        --msvc-static-runtime) STATIC_CRT="yes"; shift ;;
        --shared) SHARED="yes"; shift ;;
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
                          (sysctl -n hw.ncpu 2>/dev/null) || echo 4 )

if [ "$FORCE" = "no" ] && { [ -f "$PREFIX/lib/libTKernel.a" ] || [ -f "$PREFIX/lib/TKernel.lib" ] ||
                           [ -f "$PREFIX/lib/libTKernel.so" ] || [ -f "$PREFIX/lib/libTKernel.dylib" ]; }; then
    echo "== OpenCASCADE $OCCT_VERSION already installed in $PREFIX, skipping =="
    exit 0
fi

LIBTYPE=Static
[ "$SHARED" = "no" ] || LIBTYPE=Shared
echo "== OpenCASCADE $OCCT_VERSION ($LIBTYPE) into $PREFIX =="
occt_source "$WORK" "$WORK/occt-$OCCT_VERSION"

set --
case "$(uname -s)" in
    Darwin)
        set -- -DCMAKE_OSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-13.3}"
        [ "$UNIVERSAL" = "no" ] || set -- "$@" "-DCMAKE_OSX_ARCHITECTURES=x86_64;arm64"
        ;;
    MINGW*|MSYS*|CYGWIN*)
        # OCCT's CMake predates CMP0091, so the runtime choice needs the policy.
        if [ "$STATIC_CRT" = "yes" ]; then
            set -- -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded
        fi
        ;;
esac
occt_configure "$WORK/occt-$OCCT_VERSION" "$WORK/build" "$PREFIX" "$LIBTYPE" \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON "$@"
cmake --build "$WORK/build" --config Release -j "$JOBS"
cmake --install "$WORK/build" --config Release
occt_licenses "$WORK/occt-$OCCT_VERSION" "$PREFIX"
rm -rf "$WORK/build" "$WORK/occt-$OCCT_VERSION"
echo "== done: -DOpenCASCADE_DIR=$PREFIX/lib/cmake/opencascade =="
