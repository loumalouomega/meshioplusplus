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
#   ./configure-wasm.sh --full --build      # the lean pair and the full pair
#   ./configure-wasm.sh --full-only --build # only meshioplusplus_wasm_full[_mt]
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
#
# --full adds a second pair of artifacts, meshioplusplus_wasm_full and
# meshioplusplus_wasm_full_mt: the same core and JS API, linked with every
# library the toolchain cross-compiles -- OpenCASCADE (STEP/IGES reading and
# fitCad), zstd and lz4 (VTK XML block compression), Emscripten's bzip2
# port (libMesh .bz2), and the serial KaHIP kaffpa partitioner. OpenCASCADE
# is LGPL-2.1 with its exception, so the full pair is not plain MIT
# (doc/wasm.md#licence). Its libraries come from
# build-wasm-deps.sh --full, one prefix per threading flavour, built
# automatically when missing (OpenCASCADE takes the better part of an hour).

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
WITH_LEAN="yes"
WITH_FULL="no"
WITH_OCCT="ON"

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
  --full                         also build the full pair
                                  (meshioplusplus_wasm_full[_mt]): OpenCASCADE,
                                  zstd, lz4, bzip2 and KaHIP linked in
  --full-only                    build only the full pair
  --without-occt                 leave OpenCASCADE out of the full pair (codecs
                                  only; the artifact then stays MIT)
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
        --full) WITH_FULL="yes"; shift ;;
        --full-only) WITH_FULL="yes"; WITH_LEAN="no"; shift ;;
        --without-occt) WITH_OCCT="OFF"; shift ;;
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
VARIANTS=""
if [ "$WITH_LEAN" = "yes" ]; then
    VARIANTS="seq"
    [ "$SEQ_ONLY" = "yes" ] || VARIANTS="$VARIANTS mt"
fi
if [ "$WITH_FULL" = "yes" ]; then
    VARIANTS="$VARIANTS full"
    [ "$SEQ_ONLY" = "yes" ] || VARIANTS="$VARIANTS full-mt"
fi

# The full libraries' prefix for one threading flavour, built when missing.
full_deps_prefix() {  # full_deps_prefix <seq|mt>
    occt_arg=""
    [ "$WITH_OCCT" = "ON" ] || occt_arg="--without-occt"
    # shellcheck disable=SC2086
    fprefix=$("$SCRIPT_DIR/build-wasm-deps.sh" $occt_arg --print-full-prefix "$1")
    if [ ! -f "$fprefix/lib/libzstd.a" ] || [ ! -f "$fprefix/lib/liblz4.a" ] ||
       { [ "$WITH_OCCT" = "ON" ] && [ ! -f "$fprefix/lib/libTKernel.a" ]; }; then
        echo "== building the full artifact's libraries ($1; one-off, OpenCASCADE takes long) ==" >&2
        # shellcheck disable=SC2086
        "$SCRIPT_DIR/build-wasm-deps.sh" --full-only --threads "$1" $occt_arg >&2
    fi
    echo "$fprefix"
}

# configure_variant <seq|mt>: configure (and, with --build, build + copy) one
# variant into its own tree. The two variants MUST be separate trees: -pthread
# and the OpenMP macro are whole-translation-unit compile-time properties of the
# shared core object library, so a single tree cannot hold both.
configure_variant() {
    variant="$1"
    case "$variant" in
        mt|full-mt)
            suffix="-mt"
            backend_args="-DMESHIOPLUSPLUS_PARALLEL_BACKEND=OPENMP -DMESHIOPLUSPLUS_WASM_THREADS=ON"
            out_suffix="_mt"
            zlib_cflags="-pthread"
            flavour="mt" ;;
        *)
            suffix=""
            backend_args="-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ"
            out_suffix=""
            zlib_cflags=""
            flavour="seq" ;;
    esac
    full_args=""
    variant_root_arg="$FIND_ROOT_ARG"
    case "$variant" in
        full|full-mt)
            build_dir="$SCRIPT_DIR/wasm-${BUILD_TYPE_LC}-full${suffix}"
            out_base="meshioplusplus_wasm_full${out_suffix}"
            fprefix=$(full_deps_prefix "$flavour")
            full_args="-DMESHIOPLUSPLUS_WASM_FULL=ON -DMESHIOPLUSPLUS_WITH_OCCT=$WITH_OCCT"
            full_args="$full_args -DMESHIOPLUSPLUS_WITH_ZSTD=ON -DMESHIOPLUSPLUS_WITH_LZ4=ON"
            full_args="$full_args -DMESHIOPLUSPLUS_WITH_BZIP2=ON"
            full_args="$full_args -DMESHIOPLUSPLUS_WITH_KAHIP=ON -DKAHIP_ROOT=$fprefix"
            if [ -n "$DEPS_PREFIX" ]; then
                variant_root_arg="-DCMAKE_FIND_ROOT_PATH=$DEPS_PREFIX;$fprefix"
            else
                variant_root_arg="-DCMAKE_FIND_ROOT_PATH=$fprefix"
            fi ;;
        *)
            build_dir="$SCRIPT_DIR/wasm-${BUILD_TYPE_LC}${suffix}"
            out_base="meshioplusplus_wasm${out_suffix}" ;;
    esac

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
    case "$variant" in
        full*) echo "  full:      OCCT $WITH_OCCT, zstd, lz4, bzip2, KaHIP ($fprefix)" ;;
    esac
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
        $variant_root_arg \
        $backend_args \
        $full_args \
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
        # The full pair links the bzip2 port too, with the same race.
        case "$variant" in
            full*) EMCC_CFLAGS="$zlib_cflags" embuilder build bzip2 ;;
        esac
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
