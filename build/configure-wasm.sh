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
#   ./configure-wasm.sh --without-hdf5 --without-netcdf --build   # small build
#
# Unlike configure.sh (the native Python-extension build), this script never
# touches Python/pybind11 (-DMESHIOPLUSPLUS_BUILD_PYTHON=OFF) and always
# configures the sequential parallel backend (OpenMP/TBB/STL's parallel STL
# all have no meaningful WASM story today; SEQ is fast enough for the format
# parsers this target ships and keeps the configure deterministic). The mesh
# backend is pinned to NATIVE (-DMESHIOPLUSPLUS_MESH_BACKEND=NATIVE): the
# fastest in-memory structure — canonical Float64/Int64 storage means the
# embind boundary's typed arrays convert with no dtype dispatch, and the JS
# API shape is unchanged.
#
# HDF5 and netCDF are ON by default here, as they are for a native build, so
# the published npm artifact reads and writes CGNS/H5M/HMF/MED/Exodus and
# XDMF's Format="HDF" data path like any other build. They are not system
# libraries on this target: build-wasm-deps.sh produces a wasm32 static build
# of both into a self-contained prefix, which this script puts on
# CMAKE_FIND_ROOT_PATH (the required knob -- the Emscripten toolchain sets
# CMAKE_FIND_ROOT_PATH_MODE_PACKAGE/LIBRARY/INCLUDE to ONLY, so a plain
# CMAKE_PREFIX_PATH or HINTS is ignored). If the prefix is missing it is built
# automatically, once; that takes several minutes, so pass --deps-prefix to
# reuse one you already have, or --without-hdf5 for a fast, small build.

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_DIR=$(dirname -- "$SCRIPT_DIR")

BUILD_TYPE="Release"
WITH_ZLIB="ON"
WITH_HDF5="ON"
WITH_NETCDF="ON"
WITH_CGNSLIB="ON"
WITH_GIDPOST="ON"
DEPS_PREFIX=""
DO_BUILD="no"
SEQ_ONLY="no"

usage() {
    cat <<EOF
Usage: $0 [options]
  --build-type <type>            CMake build type (default: Release)
  --with-zlib / --without-zlib   VTU/XDMF zlib compression via Emscripten's
                                  built-in port (default: on)
  --with-hdf5 / --without-hdf5   CGNS/H5M/HMF/MED + XDMF-HDF (default: on)
  --with-netcdf / --without-netcdf
  --with-cgnslib / --without-cgnslib
                                  the CGNS MLL backend: ADF containers and the
                                  CGNS 3.x NGON_n layout. Everything else in
                                  CGNS works without it (doc/formats/cgns.md).
                                 Exodus (default: on; needs HDF5)
  --with-gidpost / --without-gidpost
                                  GiD postprocess writer, vendored gidpost
                                  (default: on; needs zlib)
  --deps-prefix <dir>            prefix holding the wasm32 HDF5/netCDF build
                                  (default: build-wasm-deps.sh --print-prefix)
  --seq-only                     build only the sequential variant (skip the
                                  threaded meshioplusplus_wasm_mt); default is
                                  to build both, since the npm package ships both
  --build                        run the build after configuring
  -h, --help                     this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --build-type) BUILD_TYPE="$2"; shift 2 ;;
        --with-zlib) WITH_ZLIB="ON"; shift ;;
        --without-zlib) WITH_ZLIB="OFF"; shift ;;
        --with-hdf5) WITH_HDF5="ON"; shift ;;
        --without-hdf5) WITH_HDF5="OFF"; WITH_NETCDF="OFF"; WITH_CGNSLIB="OFF"; shift ;;
        --with-cgnslib) WITH_CGNSLIB="ON"; shift ;;
        --without-cgnslib) WITH_CGNSLIB="OFF"; shift ;;
        --with-netcdf) WITH_NETCDF="ON"; shift ;;
        --without-netcdf) WITH_NETCDF="OFF"; shift ;;
        --with-gidpost) WITH_GIDPOST="ON"; shift ;;
        --without-gidpost) WITH_GIDPOST="OFF"; shift ;;
        --deps-prefix) DEPS_PREFIX="$2"; shift 2 ;;
        --seq-only) SEQ_ONLY="yes"; shift ;;
        --build) DO_BUILD="yes"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# netCDF-4 is layered on HDF5; asking for one without the other cannot work.
if [ "$WITH_NETCDF" = "ON" ] && [ "$WITH_HDF5" = "OFF" ]; then
    echo "error: --with-netcdf needs HDF5 (netCDF-4 is an HDF5 container)." >&2
    exit 1
fi

# The HDF5-backed writers all compress (H5Pset_deflate, gzip level 4), so an
# HDF5 build without zlib would produce files it cannot write.
if [ "$WITH_HDF5" = "ON" ] && [ "$WITH_ZLIB" = "OFF" ]; then
    echo "error: --with-hdf5 needs zlib (the MED/CGNS/H5M/HMF writers deflate)." >&2
    echo "       Use --without-hdf5 as well, or drop --without-zlib." >&2
    exit 1
fi

# gidpostInt.h includes <zlib.h> unconditionally and the binary flavour is
# always deflated, so the GiD writer cannot be built without zlib.
if [ "$WITH_GIDPOST" = "ON" ] && [ "$WITH_ZLIB" = "OFF" ]; then
    echo "error: --with-gidpost needs zlib (gidpost deflates unconditionally)." >&2
    echo "       Use --without-gidpost as well, or drop --without-zlib." >&2
    exit 1
fi

if ! command -v emcmake >/dev/null 2>&1; then
    echo "error: emcmake not found on PATH." >&2
    echo "Install the Emscripten SDK and 'source \$EMSDK/emsdk_env.sh' first:" >&2
    echo "  https://emscripten.org/docs/getting_started/downloads.html" >&2
    exit 1
fi

BUILD_TYPE_LC=$(echo "$BUILD_TYPE" | tr '[:upper:]' '[:lower:]')

GENERATOR=""
if command -v ninja >/dev/null 2>&1; then
    GENERATOR="-G Ninja"
fi

# Locate (and if necessary produce) the wasm32 HDF5/netCDF prefix. Asking the
# deps script for the default path rather than recomputing it here keeps the
# version pins in exactly one place. The prefix is shared by both variants
# (HDF5/netCDF are not affected by the pthreads/OpenMP switch).
FIND_ROOT_ARG=""
if [ "$WITH_HDF5" = "ON" ]; then
    [ -n "$DEPS_PREFIX" ] || DEPS_PREFIX=$("$SCRIPT_DIR/build-wasm-deps.sh" --print-prefix)
    if [ ! -f "$DEPS_PREFIX/lib/libhdf5.a" ] ||
       { [ "$WITH_NETCDF" = "ON" ] && [ ! -f "$DEPS_PREFIX/lib/libnetcdf.a" ]; } ||
       { [ "$WITH_CGNSLIB" = "ON" ] && [ ! -f "$DEPS_PREFIX/lib/libcgns.a" ]; }; then
        echo "== building the wasm HDF5/netCDF dependencies (one-off, several minutes) =="
        deps_args="--prefix $DEPS_PREFIX"
        [ "$WITH_NETCDF" = "ON" ] || deps_args="$deps_args --without-netcdf"
        [ "$WITH_CGNSLIB" = "ON" ] || deps_args="$deps_args --without-cgnslib"
        # shellcheck disable=SC2086
        "$SCRIPT_DIR/build-wasm-deps.sh" $deps_args
        echo
    fi
    FIND_ROOT_ARG="-DCMAKE_FIND_ROOT_PATH=$DEPS_PREFIX"
fi

# Two artifacts ship in the npm package: the sequential build
# (meshioplusplus_wasm) and the threaded OpenMP/pthreads build
# (meshioplusplus_wasm_mt). The loader (src/wasm/src/index.mjs) auto-selects
# the threaded one where the page is cross-origin isolated and falls back to the
# sequential one otherwise, so both must be present. --seq-only skips the
# threaded one for a fast, small, header-free build.
VARIANTS="seq mt"
[ "$SEQ_ONLY" = "yes" ] && VARIANTS="seq"

# configure_variant <seq|mt>: configure (and, with --build, build + copy) one
# variant into its own tree. The two variants MUST be separate trees: -pthread
# and the OpenMP macro are whole-translation-unit compile-time properties of the
# shared core object library, so a single tree cannot hold both.
configure_variant() {
    variant="$1"
    if [ "$variant" = "mt" ]; then
        build_dir="$SCRIPT_DIR/wasm-${BUILD_TYPE_LC}-mt"
        backend_args="-DMESHIOPLUSPLUS_PARALLEL_BACKEND=OPENMP -DMESHIOPLUSPLUS_WASM_THREADS=ON"
        out_base="meshioplusplus_wasm_mt"
        zlib_cflags="-pthread"
    else
        build_dir="$SCRIPT_DIR/wasm-${BUILD_TYPE_LC}"
        backend_args="-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ"
        out_base="meshioplusplus_wasm"
        zlib_cflags=""
    fi

    echo "== meshio++ WASM configure ($variant) =="
    echo "  source:    $SOURCE_DIR"
    echo "  build:     $build_dir"
    echo "  type:      $BUILD_TYPE"
    echo "  variant:   $variant ($out_base)"
    echo "  zlib:      $WITH_ZLIB"
    echo "  hdf5:      $WITH_HDF5"
    echo "  netcdf:    $WITH_NETCDF"
    echo "  cgnslib:   $WITH_CGNSLIB"
    echo "  gidpost:   $WITH_GIDPOST"
    [ "$WITH_HDF5" = "OFF" ] || echo "  deps:      $DEPS_PREFIX"
    echo "  emcc:      $(command -v emcc)"
    echo

    # CMAKE_FIND_ROOT_PATH rather than CMAKE_PREFIX_PATH: the Emscripten
    # toolchain re-roots every find_package/find_library/find_path at the
    # sysroot (CMAKE_FIND_ROOT_PATH_MODE_* = ONLY), and appends the sysroot to
    # whatever CMAKE_FIND_ROOT_PATH already holds -- so passing our prefix here
    # adds it as a second root instead of replacing the toolchain's.
    # shellcheck disable=SC2086
    emcmake cmake $GENERATOR \
        -S "$SOURCE_DIR" -B "$build_dir" \
        -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
        $FIND_ROOT_ARG \
        $backend_args \
        -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
        -DMESHIOPLUSPLUS_BUILD_WASM=ON \
        -DMESHIOPLUSPLUS_MESH_BACKEND=NATIVE \
        -DMESHIOPLUSPLUS_WITH_HDF5="$WITH_HDF5" \
        -DMESHIOPLUSPLUS_WITH_NETCDF="$WITH_NETCDF" \
        -DMESHIOPLUSPLUS_WITH_CGNSLIB="$WITH_CGNSLIB" \
        -DMESHIOPLUSPLUS_WITH_ZLIB="$WITH_ZLIB" \
        -DMESHIOPLUSPLUS_WITH_GIDPOST="$WITH_GIDPOST"

    if [ "$WITH_ZLIB" = "ON" ]; then
        echo
        echo "== warming Emscripten zlib port cache ($variant) =="
        # -sUSE_ZLIB=1 makes every translation unit trigger a build of the
        # bundled zlib port on first use. Left to a parallel `-j` build, Ninja's
        # dependency-scan step (emscan-deps) launches many em++ invocations at
        # once, and on a cold cache they race to build/lock that same port,
        # aborting with "attempt to lock the cache while a parent process is
        # holding the lock (sanity)". Building it once, single-threaded, up
        # front avoids the race entirely. The threaded variant needs the
        # -pthread port variant, which is a distinct cache entry, so it is
        # warmed with EMCC_CFLAGS=-pthread.
        EMCC_CFLAGS="$zlib_cflags" embuilder build zlib
    fi

    echo
    echo "== next steps ($variant) =="
    echo "  emmake cmake --build \"$build_dir\" -j"
    echo "  cp \"$build_dir\"/$out_base.{mjs,wasm} \"$SOURCE_DIR/src/wasm/dist/\""

    if [ "$DO_BUILD" = "yes" ]; then
        echo
        echo "== building ($variant) =="
        emmake cmake --build "$build_dir" -j
        mkdir -p "$SOURCE_DIR/src/wasm/dist"
        cp "$build_dir/$out_base.mjs" "$SOURCE_DIR/src/wasm/dist/"
        cp "$build_dir/$out_base.wasm" "$SOURCE_DIR/src/wasm/dist/"
        echo "copied $out_base.{mjs,wasm} to $SOURCE_DIR/src/wasm/dist/"
    fi
    echo
}

for v in $VARIANTS; do
    configure_variant "$v"
done

if [ "$DO_BUILD" != "yes" ]; then
    echo "== next steps (both variants) =="
    echo "  re-run with --build, then:"
    echo "  node \"$SOURCE_DIR/tests/wasm/smoke.mjs\""
    echo "  (cd \"$SOURCE_DIR/src/wasm\" && npm pack)"
fi
