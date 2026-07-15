#!/bin/sh
# configure-wasm.sh — configure (and optionally build) the meshio++
# WebAssembly target (@meshioplusplus/wasm) on Linux/macOS. Run from
# anywhere; the build tree is created next to this script
# (build/wasm-<build-type>). Requires the Emscripten SDK (emsdk) on PATH --
# see https://emscripten.org/docs/getting_started/downloads.html or:
#
#   git clone https://github.com/emscripten-core/emsdk.git
#   cd emsdk && ./emsdk install latest && ./emsdk activate latest
#   source ./emsdk_env.sh
#
#   ./configure-wasm.sh --build
#   ./configure-wasm.sh --without-zlib --build-type RelWithDebInfo --build
#
# Unlike configure.sh (the native Python-extension build), this script never
# touches Python/pybind11 (-DMESHIOPLUSPLUS_BUILD_PYTHON=OFF) and always
# configures the sequential parallel backend (OpenMP/TBB/STL's parallel STL
# all have no meaningful WASM story today; SEQ is fast enough for the format
# parsers this target ships and keeps the configure deterministic). HDF5 and
# netCDF are also off unconditionally: porting them to WASM is a separate,
# much larger undertaking (see doc/wasm.md) -- the CGNS/H5M/HMF/MED/Exodus
# formats are simply absent from bindings_js/js_bindings.cpp's dispatch
# tables regardless of these flags.

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_DIR=$(dirname -- "$SCRIPT_DIR")

BUILD_TYPE="Release"
WITH_ZLIB="ON"
DO_BUILD="no"

usage() {
    cat <<EOF
Usage: $0 [options]
  --build-type <type>          CMake build type (default: Release)
  --with-zlib / --without-zlib VTU/XDMF zlib compression via Emscripten's
                                built-in port (default: on)
  --build                      run the build after configuring
  -h, --help                   this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --build-type) BUILD_TYPE="$2"; shift 2 ;;
        --with-zlib) WITH_ZLIB="ON"; shift ;;
        --without-zlib) WITH_ZLIB="OFF"; shift ;;
        --build) DO_BUILD="yes"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if ! command -v emcmake >/dev/null 2>&1; then
    echo "error: emcmake not found on PATH." >&2
    echo "Install the Emscripten SDK and 'source \$EMSDK/emsdk_env.sh' first:" >&2
    echo "  https://emscripten.org/docs/getting_started/downloads.html" >&2
    exit 1
fi

BUILD_DIR="$SCRIPT_DIR/wasm-$(echo "$BUILD_TYPE" | tr '[:upper:]' '[:lower:]')"

GENERATOR=""
if command -v ninja >/dev/null 2>&1; then
    GENERATOR="-G Ninja"
fi

echo "== meshio++ WASM configure =="
echo "  source:    $SOURCE_DIR"
echo "  build:     $BUILD_DIR"
echo "  type:      $BUILD_TYPE"
echo "  zlib:      $WITH_ZLIB"
echo "  emcc:      $(command -v emcc)"
echo

# shellcheck disable=SC2086
emcmake cmake $GENERATOR \
    -S "$SOURCE_DIR" -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
    -DMESHIOPLUSPLUS_BUILD_WASM=ON \
    -DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ \
    -DMESHIOPLUSPLUS_WITH_HDF5=OFF \
    -DMESHIOPLUSPLUS_WITH_NETCDF=OFF \
    -DMESHIOPLUSPLUS_WITH_ZLIB="$WITH_ZLIB"

if [ "$WITH_ZLIB" = "ON" ]; then
    echo
    echo "== warming Emscripten zlib port cache =="
    # -sUSE_ZLIB=1 makes every translation unit trigger a build of the
    # bundled zlib port on first use. Left to a parallel `-j` build, Ninja's
    # dependency-scan step (emscan-deps) launches many em++ invocations at
    # once, and on a cold cache they race to build/lock that same port,
    # aborting with "attempt to lock the cache while a parent process is
    # holding the lock (sanity)". Building it once, single-threaded, up
    # front avoids the race entirely.
    embuilder build zlib
fi

echo
echo "== next steps =="
echo "  emmake cmake --build \"$BUILD_DIR\" -j"
echo "  cp \"$BUILD_DIR\"/meshioplusplus_wasm.{mjs,wasm} \"$SOURCE_DIR/wasm/dist/\""
echo "  node \"$SOURCE_DIR/wasm/test/smoke.mjs\""
echo "  (cd \"$SOURCE_DIR/wasm\" && npm pack)"

if [ "$DO_BUILD" = "yes" ]; then
    echo
    echo "== building =="
    emmake cmake --build "$BUILD_DIR" -j
    mkdir -p "$SOURCE_DIR/wasm/dist"
    cp "$BUILD_DIR"/meshioplusplus_wasm.mjs "$SOURCE_DIR/wasm/dist/"
    cp "$BUILD_DIR"/meshioplusplus_wasm.wasm "$SOURCE_DIR/wasm/dist/"
    echo "copied build artifacts to $SOURCE_DIR/wasm/dist/"
fi
