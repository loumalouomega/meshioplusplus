#!/bin/sh
# build-wasm-deps.sh — build the wasm32-emscripten static libraries that the
# HDF5- and netCDF-backed formats (CGNS, H5M, HMF, MED, XDMF's Format="HDF",
# and Exodus) need, and install them into a self-contained prefix that
# build/configure-wasm.sh then points CMake at.
#
#   ./build-wasm-deps.sh                 # build both into the default prefix
#   ./build-wasm-deps.sh --print-prefix  # just print that prefix and exit
#   ./build-wasm-deps.sh --without-netcdf
#   ./build-wasm-deps.sh --full          # also the full artifact's libraries
#   ./build-wasm-deps.sh --full-only --threads seq
#   ./build-wasm-deps.sh --print-full-prefix mt
#
# The full artifact (meshioplusplus_wasm_full, configure-wasm.sh --full) links
# five more libraries: OpenCASCADE (STEP/IGES reading and fit_cad), zstd and
# lz4 (VTK XML block compression), the serial KaHIP kaffpa partitioner and
# Netgen's nglib (the tetrahedralize backend that keeps the surface).
# bzip2 is Emscripten's own port, like zlib, and needs nothing here. They
# install into their own prefix, one per threading flavour: OCCT's reference
# counts are std::atomic, and an object compiled without -pthread has its
# atomics lowered to plain loads and stores and is refused by the threaded
# (--shared-memory) link. The lean prefix and its CI cache are untouched by
# --full.
#
# Requires the Emscripten SDK (emsdk) on PATH -- the same one used for the
# meshio++ build itself, since a static archive is only guaranteed to link
# against objects produced by a compatible toolchain. See
# https://emscripten.org/docs/getting_started/downloads.html or:
#
#   git clone https://github.com/emscripten-core/emsdk.git
#   cd emsdk && ./emsdk install latest && ./emsdk activate latest
#   source ./emsdk_env.sh
#
# Why a shell script and not CMake: meshio++'s CMake never downloads anything
# (the KaHIP rule) -- it only ever *finds* an already-built dependency. This
# script is the documented, reproducible way to produce one for the wasm
# target, mirroring how CI source-builds KaHIP into a cached prefix.
#
# Both libraries are built with -Oz rather than -O3: they are I/O- and
# metadata-bound (B-tree walks, dataspace bookkeeping), not hot numeric loops,
# so the size saving is close to free on a target where every byte is shipped
# over the network.
#
# The prefix is NOT relocatable -- netCDF's installed CMake config records the
# absolute path of the HDF5 libraries it was configured against. Move it and
# re-run this script rather than moving the directory.

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

# Pinned upstream versions. HDF5 >= 1.14.0 is required, not merely preferred:
# it is the release that removed H5detect/H5make_libsettings (their work now
# happens at library startup and via template files), which is exactly what
# used to make cross-compiling HDF5 painful.
HDF5_VERSION="1.14.6"
NETCDF_VERSION="4.9.3"
CGNS_VERSION="4.5.2"
HDF5_SHA256="e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b"
NETCDF_SHA256="990f46d49525d6ab5dc4249f8684c6deeaf54de6fec63a187e9fb382cc0ffdff"
CGNS_SHA256="95075e1fd0b51d97b1b96b73ebe03b1a551fbcc9cd2b2b6f487ccccedcff5964"
# The commits the netCDF and CGNS tags peel to. GitHub serves a tag's archive
# from codeload.github.com, which some proxies refuse; the fetch then falls
# back to a shallow git clone of the tag, verified against this commit.
NETCDF_COMMIT="3a6b1debf1557f07b606ce3653e44f0d711203be"
CGNS_COMMIT="3a972d92559787a49c5f9a54974ea1946ecb9d21"

# The full artifact's libraries. OpenCASCADE's pin, source and toolkit list are
# build/occt-recipe.sh's, shared with the native static build the wheels and
# CLI binaries link (build/build-occt.sh). zstd and lz4 are their release
# tarballs, whose checksums match the .sha256 files published beside them.
# KaHIP is its release tarball, checksummed the same way (or, where the archive
# is refused, the clone of its tag verified against KAHIP_COMMIT); the serial
# kaffpa interface only (PARHIP off, no MPI, no OpenMP under Emscripten),
# matching the native KaHIP integration's link surface. Netgen is a shallow
# clone of its tag, refused unless it is NETGEN_COMMIT, as OpenCASCADE is: its
# nglib and ngcore only (no GUI, Python, OCC, MPI or CGNS), the tetrahedralize
# backend, and the same 6.2.2401 the native CI job links (Ubuntu 24.04's
# libnglib-dev).
OCCT_RECIPE_DIR="$SCRIPT_DIR"
. "$SCRIPT_DIR/occt-recipe.sh"
ZSTD_VERSION="1.5.7"
ZSTD_SHA256="eb33e51f49a15e023950cd7825ca74a4a2b43db8354825ac24fc1b7ee09e6fa3"
LZ4_VERSION="1.10.0"
LZ4_SHA256="537512904744b35e232912055ccf8ec66d768639ff3abe5788d90d792ec5f48b"
KAHIP_VERSION="3.25"
KAHIP_SHA256="3abad20158887b585b4d4792c35fa9023b7fd634b2bc4494334a3951b69f4842"
KAHIP_COMMIT="2e86d3da2065bcf6ddbf6fb84d2b52f33db1061f"
NETGEN_VERSION="6.2.2401"
NETGEN_COMMIT="2ff886457fe899685d9d826c8a29f58a31889dda"

WITH_HDF5="yes"
WITH_NETCDF="yes"
WITH_CGNS="yes"
PREFIX=""
PRINT_PREFIX="no"
FORCE="no"
JOBS=""
WITH_LEAN="yes"
WITH_FULL="no"
WITH_OCCT="yes"
THREADS="both"
PRINT_FULL_PREFIX=""

usage() {
    cat <<EOF
Usage: $0 [options]
  --prefix <dir>            install prefix (default: see --print-prefix)
  --print-prefix            print the default prefix and exit
  --without-hdf5            skip HDF5 (implies --without-netcdf)
  --without-netcdf          skip netCDF
  --without-cgnslib         skip cgnslib (the CGNS MLL backend)
  --full                    also build the full artifact's libraries (OCCT,
                            zstd, lz4, KaHIP) into their own prefix(es)
  --full-only               build only those (skip HDF5/netCDF/cgnslib)
  --threads <seq|mt|both>   which flavours of the full libraries (default: both)
  --without-occt            leave OpenCASCADE out of the full prefix
  --print-full-prefix <seq|mt>
                            print the full libraries' prefix and exit
  --force                   rebuild even if the prefix already looks complete
  --jobs <n>                parallel build jobs (default: all cores)
  -h, --help                this help

Pinned versions: HDF5 $HDF5_VERSION, netcdf-c $NETCDF_VERSION, cgnslib $CGNS_VERSION,
                 OCCT $OCCT_VERSION, zstd $ZSTD_VERSION, lz4 $LZ4_VERSION, KaHIP $KAHIP_VERSION (full only)
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) PREFIX="$2"; shift 2 ;;
        --print-prefix) PRINT_PREFIX="yes"; shift ;;
        --without-hdf5) WITH_HDF5="no"; WITH_NETCDF="no"; WITH_CGNS="no"; shift ;;
        --without-cgnslib) WITH_CGNS="no"; shift ;;
        --without-netcdf) WITH_NETCDF="no"; shift ;;
        --full) WITH_FULL="yes"; shift ;;
        --full-only) WITH_FULL="yes"; WITH_LEAN="no"; shift ;;
        --threads) THREADS="$2"; shift 2 ;;
        --without-occt) WITH_OCCT="no"; shift ;;
        --print-full-prefix) PRINT_FULL_PREFIX="$2"; shift 2 ;;
        --force) FORCE="yes"; shift ;;
        --jobs) JOBS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# The prefix name carries both versions so a version bump lands in a fresh
# directory instead of half-overwriting the old one (and so CI's cache key can
# simply be the directory name).
DEPS_ROOT="$SCRIPT_DIR/wasm-deps"
[ -n "$PREFIX" ] || PREFIX="$DEPS_ROOT/hdf5-$HDF5_VERSION-netcdf-$NETCDF_VERSION-cgns-$CGNS_VERSION"

# The full libraries' prefix, per flavour. Its name carries every pinned
# version for the same reason the lean one does; --without-occt gets its own.
full_prefix() {  # full_prefix <seq|mt>
    name="full"
    [ "$WITH_OCCT" = "yes" ] && name="$name-occt-$OCCT_VERSION"
    name="$name-zstd-$ZSTD_VERSION-lz4-$LZ4_VERSION-kahip-$KAHIP_VERSION-netgen-$NETGEN_VERSION"
    [ "$1" = "mt" ] && name="$name-mt"
    echo "$DEPS_ROOT/$name"
}

if [ "$PRINT_PREFIX" = "yes" ]; then
    echo "$PREFIX"
    exit 0
fi
case "$PRINT_FULL_PREFIX" in
    "") ;;
    seq|mt) full_prefix "$PRINT_FULL_PREFIX"; exit 0 ;;
    *) echo "error: --print-full-prefix takes seq or mt" >&2; exit 1 ;;
esac
case "$THREADS" in
    seq) FULL_FLAVOURS="seq" ;;
    mt) FULL_FLAVOURS="mt" ;;
    both) FULL_FLAVOURS="seq mt" ;;
    *) echo "error: --threads takes seq, mt or both" >&2; exit 1 ;;
esac

if [ "$WITH_LEAN" = "no" ]; then
    WITH_HDF5="no"; WITH_NETCDF="no"; WITH_CGNS="no"
fi
if [ "$WITH_HDF5" = "no" ] && [ "$WITH_NETCDF" = "no" ] && [ "$WITH_FULL" = "no" ]; then
    echo "nothing to do (--without-hdf5 --without-netcdf)"
    exit 0
fi

if ! command -v emcmake >/dev/null 2>&1; then
    echo "error: emcmake not found on PATH." >&2
    echo "Install the Emscripten SDK and 'source \$EMSDK/emsdk_env.sh' first:" >&2
    echo "  https://emscripten.org/docs/getting_started/downloads.html" >&2
    exit 1
fi

for tool in curl tar patch cmake git; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "error: '$tool' not found on PATH (needed to fetch, unpack, patch and build)." >&2
        exit 1
    }
done

[ -n "$JOBS" ] || JOBS=$( (command -v nproc >/dev/null 2>&1 && nproc) || echo 4 )

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | cut -d' ' -f1
    else
        shasum -a 256 "$1" | cut -d' ' -f1  # macOS
    fi
}

# Download into a version-independent cache next to the prefixes, so bumping a
# version (or --force) never re-downloads what is already there.
DL_DIR="$DEPS_ROOT/_dl"
WORK_DIR="$DEPS_ROOT/_work"
mkdir -p "$DL_DIR" "$WORK_DIR"

fetch() {  # fetch <url> <file> <sha256>; returns 1 when the download fails
    url="$1"; file="$DL_DIR/$2"; want="$3"
    if [ -f "$file" ] && [ "$(sha256_of "$file")" = "$want" ]; then
        echo "  cached: $2"
        return 0
    fi
    echo "  downloading: $2"
    if ! curl -fsSL -o "$file.part" "$url"; then
        rm -f "$file.part"
        echo "  download failed: $url" >&2
        return 1
    fi
    got=$(sha256_of "$file.part")
    if [ "$got" != "$want" ]; then
        rm -f "$file.part"
        echo "error: checksum mismatch for $2" >&2
        echo "  expected $want" >&2
        echo "  got      $got" >&2
        exit 1
    fi
    mv "$file.part" "$file"
}

# fetch_git <url> <tag> <commit> <dir>: a shallow clone of <tag> into
# $DL_DIR/<dir>, refused unless it checks out exactly <commit>; a cached clone
# at that commit is reused.
fetch_git() {
    url="$1"; tag="$2"; commit="$3"; dir="$DL_DIR/$4"
    if [ -d "$dir/.git" ] && [ "$(git -C "$dir" rev-parse HEAD)" = "$commit" ]; then
        echo "  cached: $4 ($tag)"
        return 0
    fi
    echo "  cloning: $url $tag"
    rm -rf "$dir"
    git -c advice.detachedHead=false clone --quiet --depth 1 --branch "$tag" "$url" "$dir"
    got=$(git -C "$dir" rev-parse HEAD)
    if [ "$got" != "$commit" ]; then
        rm -rf "$dir"
        echo "error: $url $tag is $got, expected $commit" >&2
        exit 1
    fi
}

# unpack_tag <archive-url> <file> <sha256> <git-url> <tag> <commit> <srcdir>:
# the checksummed tag archive unpacked into $WORK_DIR, or, when the archive
# cannot be downloaded, the commit-verified clone of the same tag copied there
# (the two hold the same tree).
unpack_tag() {
    if fetch "$1" "$2" "$3"; then
        tar -xzf "$DL_DIR/$2" -C "$WORK_DIR"
    else
        echo "  archive unavailable, falling back to git"
        fetch_git "$4" "$5" "$6" "$7.git"
        rm -rf "$WORK_DIR/$7"
        mkdir -p "$WORK_DIR/$7"
        (cd "$DL_DIR/$7.git" && git archive --format=tar HEAD) | tar -xf - -C "$WORK_DIR/$7"
    fi
}

# Patches live in build/patches/<name>-*.patch and are applied in shell glob
# order. They are vendored rather than fetched so the build is reproducible
# offline and the diff is reviewable in-tree; each one carries its rationale
# in the patched comment itself.
apply_patches() {  # apply_patches <name> <srcdir>
    for patchfile in "$SCRIPT_DIR/patches/$1-"*.patch; do
        [ -e "$patchfile" ] || continue
        echo "  applying $(basename "$patchfile")"
        (cd "$2" && patch -p1 --forward --silent < "$patchfile")
    done
}

SYSROOT="$(em-config CACHE)/sysroot"

echo "== meshio++ WASM dependency build =="
echo "  prefix:   $PREFIX"
echo "  hdf5:     $([ "$WITH_HDF5" = yes ] && echo "$HDF5_VERSION" || echo "(skipped)")"
echo "  netcdf:   $([ "$WITH_NETCDF" = yes ] && echo "$NETCDF_VERSION" || echo "(skipped)")"
echo "  cgnslib:  $([ "$WITH_CGNS" = yes ] && echo "$CGNS_VERSION" || echo "(skipped)")"
if [ "$WITH_FULL" = "yes" ]; then
    echo "  full:     zstd $ZSTD_VERSION, lz4 $LZ4_VERSION, KaHIP $KAHIP_VERSION$([ "$WITH_OCCT" = yes ] && echo ", OCCT $OCCT_VERSION") ($FULL_FLAVOURS)"
fi
echo "  emcc:     $(command -v emcc)"
echo "  sysroot:  $SYSROOT"
echo "  jobs:     $JOBS"
echo

# Emscripten ships zlib as a *port*: -sUSE_ZLIB=1 supplies both <zlib.h> and
# the implementation, but only once the port has been materialised into the
# sysroot. Build it up front, single-threaded, for the same reason
# configure-wasm.sh does: a parallel build racing to populate a cold port
# cache aborts with "attempt to lock the cache while a parent process is
# holding the lock".
echo "== warming Emscripten zlib port cache =="
embuilder build zlib
ZLIB_INCLUDE_DIR="$SYSROOT/include"
ZLIB_LIBRARY="$SYSROOT/lib/wasm32-emscripten/libz.a"
if [ ! -f "$ZLIB_LIBRARY" ]; then
    echo "error: expected the zlib port at $ZLIB_LIBRARY after 'embuilder build zlib'" >&2
    exit 1
fi

# HDF5's and netCDF's CMake both run try_compile/try_run configure probes;
# under emcmake those work because the toolchain sets
# CMAKE_CROSSCOMPILING_EMULATOR to node.
export CFLAGS="${CFLAGS:-} -Oz -sUSE_ZLIB=1"
export CXXFLAGS="${CXXFLAGS:-} -Oz -sUSE_ZLIB=1"

# --------------------------------------------------------------------------
# HDF5
# --------------------------------------------------------------------------
if [ "$WITH_HDF5" = "yes" ]; then
    if [ "$FORCE" = "no" ] && [ -f "$PREFIX/lib/libhdf5.a" ]; then
        echo "== HDF5 $HDF5_VERSION already installed, skipping (use --force to rebuild) =="
    else
        echo "== fetching HDF5 $HDF5_VERSION =="
        fetch "https://github.com/HDFGroup/hdf5/releases/download/hdf5_$HDF5_VERSION/hdf5-$HDF5_VERSION.tar.gz" \
              "hdf5-$HDF5_VERSION.tar.gz" "$HDF5_SHA256" || exit 1
        rm -rf "$WORK_DIR/hdf5-$HDF5_VERSION" "$WORK_DIR/build-hdf5"
        tar -xzf "$DL_DIR/hdf5-$HDF5_VERSION.tar.gz" -C "$WORK_DIR"
        apply_patches "hdf5-$HDF5_VERSION" "$WORK_DIR/hdf5-$HDF5_VERSION"

        echo "== configuring HDF5 =="
        # Everything that needs an OS service wasm does not have is off:
        # threadsafe (pthreads), parallel (MPI), the plugin loader (dlopen),
        # the ROS3/direct/mirror/subfiling VFDs (curl, O_DIRECT, sockets, MPI).
        # HL *is* built: meshio++ itself only uses the core C API, but
        # netCDF-4's dimension-scale layer needs H5DS*, and building it here is
        # far cheaper than discovering that during the netCDF link.
        emcmake cmake -S "$WORK_DIR/hdf5-$HDF5_VERSION" -B "$WORK_DIR/build-hdf5" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="$PREFIX" \
            -DBUILD_SHARED_LIBS=OFF \
            -DBUILD_STATIC_LIBS=ON \
            -DONLY_SHARED_LIBS=OFF \
            -DBUILD_TESTING=OFF \
            -DHDF5_BUILD_TOOLS=OFF \
            -DHDF5_BUILD_UTILS=OFF \
            -DHDF5_BUILD_EXAMPLES=OFF \
            -DHDF5_BUILD_CPP_LIB=OFF \
            -DHDF5_BUILD_FORTRAN=OFF \
            -DHDF5_BUILD_JAVA=OFF \
            -DHDF5_BUILD_HL_LIB=ON \
            -DHDF5_ENABLE_DEPRECATED_SYMBOLS=ON \
            -DHDF5_ENABLE_THREADSAFE=OFF \
            -DHDF5_ENABLE_PARALLEL=OFF \
            -DHDF5_ENABLE_PLUGIN_SUPPORT=OFF \
            -DHDF5_ENABLE_ROS3_VFD=OFF \
            -DHDF5_ENABLE_DIRECT_VFD=OFF \
            -DHDF5_ENABLE_MIRROR_VFD=OFF \
            -DHDF5_ENABLE_SUBFILING_VFD=OFF \
            -DHDF5_ENABLE_SZIP_SUPPORT=OFF \
            -DHDF5_ENABLE_Z_LIB_SUPPORT=ON \
            -DZLIB_INCLUDE_DIR="$ZLIB_INCLUDE_DIR" \
            -DZLIB_LIBRARY="$ZLIB_LIBRARY" \
            -DZLIB_USE_EXTERNAL=OFF \
            -DHDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16=OFF

        echo "== building HDF5 =="
        cmake --build "$WORK_DIR/build-hdf5" -j "$JOBS"
        cmake --install "$WORK_DIR/build-hdf5"

        # HDF5 bakes the literal imported-target name ZLIB::ZLIB into the
        # exported link interface of hdf5-static (CMakeFilters.cmake), but its
        # installed config package never re-creates that target. A consumer
        # therefore inherits a dangling name -- and CMake's try_compile scratch
        # projects, which netcdf-c's configure uses heavily, fail outright on
        # it ("the target was not found"). Rewrite it to the archive's absolute
        # path so the installed package is self-contained. The archive is
        # Emscripten's own zlib port, the same one our final link pulls in via
        # -sUSE_ZLIB=1, so both spellings resolve to the same bytes.
        echo "== fixing up the exported HDF5 link interface =="
        sed -i.bak "s|ZLIB::ZLIB|$ZLIB_LIBRARY|g" "$PREFIX/cmake/hdf5-targets.cmake"
        rm -f "$PREFIX/cmake/hdf5-targets.cmake.bak"
    fi
fi

# --------------------------------------------------------------------------
# netCDF-C
# --------------------------------------------------------------------------
if [ "$WITH_NETCDF" = "yes" ]; then
    if [ "$FORCE" = "no" ] && [ -f "$PREFIX/lib/libnetcdf.a" ]; then
        echo "== netCDF $NETCDF_VERSION already installed, skipping (use --force to rebuild) =="
    else
        echo "== fetching netcdf-c $NETCDF_VERSION =="
        rm -rf "$WORK_DIR/netcdf-c-$NETCDF_VERSION" "$WORK_DIR/build-netcdf"
        unpack_tag "https://github.com/Unidata/netcdf-c/archive/refs/tags/v$NETCDF_VERSION.tar.gz" \
                   "netcdf-c-$NETCDF_VERSION.tar.gz" "$NETCDF_SHA256" \
                   "https://github.com/Unidata/netcdf-c.git" "v$NETCDF_VERSION" "$NETCDF_COMMIT" \
                   "netcdf-c-$NETCDF_VERSION"
        apply_patches "netcdf-c-$NETCDF_VERSION" "$WORK_DIR/netcdf-c-$NETCDF_VERSION"

        echo "== configuring netcdf-c =="
        # Every option is passed under both spellings: 4.9.3 renamed them to a
        # NETCDF_ENABLE_* prefix (keeping the old names as deprecated aliases),
        # and an unknown -D cache entry is only a warning -- so one command
        # line works across the rename in either direction.
        #
        # What is switched off, and why it matters here: DAP/DAP4/byterange/S3
        # need libcurl, NCZarr needs curl+zip, LIBXML2=OFF falls back to the
        # bundled ezxml, and PLUGINS=OFF removes the dlopen-based filter
        # loader. MMAP=OFF drops the NC_DISKLESS mmap backend, which nothing
        # here uses and whose mmapio.c declares mremap with a 4-argument
        # signature -- against Emscripten's 5-argument libc one, that is a
        # wasm-ld "function signature mismatch" warning on every link.
        # What remains is libsrc (classic netCDF-3) plus libhdf5
        # (netCDF-4) -- exactly what src/cpp/src/formats/exodus.cpp needs,
        # since it creates with NC_CLOBBER | NC_NETCDF4 and must still read
        # classic files.
        emcmake cmake -S "$WORK_DIR/netcdf-c-$NETCDF_VERSION" -B "$WORK_DIR/build-netcdf" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="$PREFIX" \
            -DCMAKE_PREFIX_PATH="$PREFIX" \
            -DCMAKE_FIND_ROOT_PATH="$PREFIX" \
            -DBUILD_SHARED_LIBS=OFF \
            -DBUILD_TESTING=OFF \
            -DHDF5_USE_STATIC_LIBRARIES=ON \
            -DENABLE_HDF5=ON            -DNETCDF_ENABLE_HDF5=ON \
            -DENABLE_DAP=OFF            -DNETCDF_ENABLE_DAP=OFF \
            -DENABLE_DAP4=OFF           -DNETCDF_ENABLE_DAP4=OFF \
            -DENABLE_BYTERANGE=OFF      -DNETCDF_ENABLE_BYTERANGE=OFF \
            -DENABLE_NCZARR=OFF         -DNETCDF_ENABLE_NCZARR=OFF \
            -DENABLE_S3=OFF             -DNETCDF_ENABLE_S3=OFF \
            -DENABLE_LIBXML2=OFF        -DNETCDF_ENABLE_LIBXML2=OFF \
            -DENABLE_MMAP=OFF           -DNETCDF_ENABLE_MMAP=OFF \
            -DENABLE_PLUGINS=OFF        -DNETCDF_ENABLE_PLUGINS=OFF \
            -DENABLE_FILTER_TESTING=OFF -DNETCDF_ENABLE_FILTER_TESTING=OFF \
            -DENABLE_TESTS=OFF          -DNETCDF_ENABLE_TESTS=OFF \
            -DENABLE_EXAMPLES=OFF       -DNETCDF_ENABLE_EXAMPLES=OFF \
            -DBUILD_UTILITIES=OFF       -DNETCDF_BUILD_UTILITIES=OFF

        echo "== building netcdf-c =="
        cmake --build "$WORK_DIR/build-netcdf" -j "$JOBS"
        cmake --install "$WORK_DIR/build-netcdf"

        # Same class of problem as the HDF5 fixup above, one level up:
        # netCDF::netcdf's exported link interface names the imported targets
        # HDF5::HDF5 and hdf5::hdf5_hl. The first comes from module-mode
        # FindHDF5 and the second only from HDF5's own config package, so no
        # single find_package call in a consumer defines both. Replace them
        # with archive paths -- and while doing so, put hdf5_hl *before* hdf5,
        # since wasm-ld resolves archives in command-line order and the
        # high-level library depends on the core one (the exported order has
        # them the other way round, plus a duplicate).
        echo "== fixing up the exported netCDF link interface =="
        sed -i.bak \
            -e "s|INTERFACE_LINK_LIBRARIES \".*\"|INTERFACE_LINK_LIBRARIES \"$PREFIX/lib/libhdf5_hl.a;$PREFIX/lib/libhdf5.a;$ZLIB_LIBRARY\"|" \
            "$PREFIX/lib/cmake/netCDF/netCDFTargets.cmake"
        rm -f "$PREFIX/lib/cmake/netCDF/netCDFTargets.cmake.bak"
    fi
fi

# --------------------------------------------------------------------------
# cgnslib (the CGNS Mid-Level Library)
#
# meshio++ reads and writes CGNS itself, over raw HDF5, so this is strictly
# an ADDITION: it buys ADF-backed containers (which are not HDF5 at all and
# so unreachable from the hand-rolled path by construction) and the CGNS 3.x
# NGON_n/NFACE_n section layout. The 4.0 layout, and everything else, works
# without it -- see doc/formats/cgns.md.
#
# Fortran, the tools and the tests are all off: only the C MLL is linked, and
# the tools want X11/Tcl, which do not exist on this target.
# --------------------------------------------------------------------------
if [ "$WITH_CGNS" = "yes" ]; then
    if [ "$FORCE" = "no" ] && [ -f "$PREFIX/lib/libcgns.a" ]; then
        echo "== cgnslib $CGNS_VERSION already installed, skipping (use --force to rebuild) =="
    else
        echo "== fetching cgnslib $CGNS_VERSION =="
        rm -rf "$WORK_DIR/CGNS-$CGNS_VERSION" "$WORK_DIR/build-cgns"
        unpack_tag "https://github.com/CGNS/CGNS/archive/refs/tags/v$CGNS_VERSION.tar.gz" \
                   "cgns-$CGNS_VERSION.tar.gz" "$CGNS_SHA256" \
                   "https://github.com/CGNS/CGNS.git" "v$CGNS_VERSION" "$CGNS_COMMIT" \
                   "CGNS-$CGNS_VERSION"
        apply_patches "CGNS-$CGNS_VERSION" "$WORK_DIR/CGNS-$CGNS_VERSION"

        echo "== configuring cgnslib =="
        # cgsize_t is 32-bit here and there is no way to change that:
        # CGNS's own CMakeLists FORCES CGNS_ENABLE_64BIT off whenever
        # CMAKE_SIZEOF_VOID_P <= 4, which wasm32 is. That is fine -- 
        # `cgns_mll.cpp` is written in terms of `cgsize_t` throughout rather
        # than a fixed width, so it compiles and reads correctly either way,
        # and the ~2^31 element ceiling is far beyond anything a browser tab
        # loads. Do NOT "fix" this by passing CGNS_ENABLE_64BIT=ON; it is
        # accepted and then silently overridden, which is worse than not
        # passing it because it reads like a guarantee.
        # NOTE the source dir is the TOP LEVEL, not src/: the root
        # CMakeLists runs the CHECK_TYPE_SIZE calls that src/CMakeLists.txt
        # then reads to pick cgsize_t's underlying type. Configuring src/
        # directly leaves them empty and fails with "Can't find suitable
        # int64_t", which reads like a cross-compilation problem but is not.
        emcmake cmake -S "$WORK_DIR/CGNS-$CGNS_VERSION" -B "$WORK_DIR/build-cgns" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="$PREFIX" \
            -DCMAKE_PREFIX_PATH="$PREFIX" \
            -DCMAKE_FIND_ROOT_PATH="$PREFIX" \
            -DBUILD_SHARED_LIBS=OFF \
            -DCGNS_BUILD_SHARED=OFF \
            -DCGNS_USE_SHARED=OFF \
            -DCGNS_ENABLE_HDF5=ON \
            -DCGNS_ENABLE_FORTRAN=OFF \
            -DCGNS_ENABLE_PARALLEL=OFF \
            -DCGNS_ENABLE_TESTS=OFF \
            -DCGNS_BUILD_TESTING=OFF \
            -DCGNS_BUILD_CGNSTOOLS=OFF \
            -DHDF5_NEED_ZLIB=ON \
            -DHDF5_INCLUDE_PATH="$PREFIX/include" \
            -DHDF5_LIBRARY="$PREFIX/lib/libhdf5.a"

        echo "== building cgnslib =="
        cmake --build "$WORK_DIR/build-cgns" -j "$JOBS"
        cmake --install "$WORK_DIR/build-cgns"

        # Same class of fixup as HDF5's and netCDF's above: cgnslib's exported
        # link interface names imported targets its own config package does
        # not create on this target. Rewrite them to archive paths so the
        # installed package is self-contained, hdf5 before zlib since wasm-ld
        # resolves archives in command-line order.
        echo "== fixing up the exported cgnslib link interface =="
        for f in "$PREFIX/lib/cmake/cgns"/*.cmake "$PREFIX/cmake"/cgns*.cmake; do
            [ -f "$f" ] || continue
            sed -i.bak \
                -e "s|hdf5-static|$PREFIX/lib/libhdf5.a|g" \
                -e "s|hdf5::hdf5-static|$PREFIX/lib/libhdf5.a|g" \
                -e "s|HDF5::HDF5|$PREFIX/lib/libhdf5.a|g" \
                -e "s|ZLIB::ZLIB|$ZLIB_LIBRARY|g" "$f"
            rm -f "$f.bak"
        done
    fi
fi

# $ZLIB_LIBRARY is Emscripten's zlib PORT archive, which lives under
# `em-config CACHE` -- i.e. under $EMSDK, which on CI is a fresh directory
# named after the *current job's* runner-assigned temp dir every single run.
# $PREFIX itself (the fixups above) stays valid across runs -- the checkout
# always lands at the same workspace path -- but a $ZLIB_LIBRARY baked in by
# an EARLIER run and then reused via CI's `actions/cache` on wasm-deps/ points
# at that earlier run's now-nonexistent $EMSDK, so ninja fails with "missing
# and no known rule to make it" on a cache hit. The two fixups above already
# make this run's paths correct when they actually execute (a fresh build),
# but they are skipped entirely on the "already installed" fast path above --
# so unconditionally re-point any stale libz.a reference at the CURRENT run's
# $ZLIB_LIBRARY, whether or not this run just built it. A same-run rebuild
# rewrites the already-correct path to itself, which is a no-op.
#
# The path-prefix wildcard is a POSITIVE class of characters legal in a Unix
# path (letters/digits/`_.-/`) -- not a negated "exclude quote/semicolon/
# whitespace" class. HDF5's exported static-lib interface wraps a private
# dependency like zlib in a $<LINK_ONLY:...> generator expression, and `$`/
# `<`/`:` are not special to a *negated* class, so it greedily swallows the
# `$<LINK_ONLY:` prefix too -- replacing the whole span with just
# $ZLIB_LIBRARY and stranding the expression's closing `>` right after the
# new path, which is exactly the "libz.a>" ninja reads and cannot find. The
# positive class can never start a match inside `$<LINK_ONLY:`, so the
# wrapper (and its `>`) is left alone and the generator expression stays
# well-formed.
for target_file in "$PREFIX/cmake/hdf5-targets.cmake" \
                    "$PREFIX/lib/cmake/netCDF/netCDFTargets.cmake" \
                    "$PREFIX/lib/cmake/cgns"/*.cmake \
                    "$PREFIX/cmake"/cgns*.cmake; do
    [ -f "$target_file" ] || continue
    sed -i.bak \
        "s#[A-Za-z0-9_./-]*/sysroot/lib/wasm32-emscripten/libz\\.a#${ZLIB_LIBRARY}#g" \
        "$target_file"
    rm -f "$target_file.bak"
done

# --------------------------------------------------------------------------
# The full artifact's libraries: OpenCASCADE, zstd, lz4, KaHIP, Netgen (--full)
#
# One prefix per threading flavour (see the header). Each is built with
# -fwasm-exceptions, the exception model the core compiles with: fit_cad and
# the STEP/IGES readers catch Standard_Failure across the library boundary,
# and wasm-ld refuses to mix exception models. The lean CFLAGS above (-Oz,
# the zlib port) do not apply: only Netgen takes the zlib port (its file I/O,
# which tetrahedralize does not call), and OCCT and Netgen are compute-bound,
# so they keep their Release optimisation.
# --------------------------------------------------------------------------
build_full_flavour() {  # build_full_flavour <seq|mt>
    flavour="$1"
    fprefix=$(full_prefix "$flavour")
    pt=""
    [ "$flavour" = "mt" ] && pt="-pthread"
    echo "== full libraries ($flavour) into $fprefix =="

    if [ "$FORCE" = "no" ] && [ -f "$fprefix/lib/libzstd.a" ]; then
        echo "== zstd $ZSTD_VERSION ($flavour) already installed, skipping =="
    else
        echo "== zstd $ZSTD_VERSION ($flavour) =="
        fetch "https://github.com/facebook/zstd/releases/download/v$ZSTD_VERSION/zstd-$ZSTD_VERSION.tar.gz" \
              "zstd-$ZSTD_VERSION.tar.gz" "$ZSTD_SHA256" || exit 1
        rm -rf "$WORK_DIR/zstd-$ZSTD_VERSION" "$WORK_DIR/build-zstd-$flavour"
        tar -xzf "$DL_DIR/zstd-$ZSTD_VERSION.tar.gz" -C "$WORK_DIR"
        # The decoder and the one-shot encoder only: no programs, no worker
        # threads (the core compresses block by block), no legacy formats.
        CFLAGS="-Oz $pt" emcmake cmake -S "$WORK_DIR/zstd-$ZSTD_VERSION/build/cmake" \
            -B "$WORK_DIR/build-zstd-$flavour" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="$fprefix" \
            -DZSTD_BUILD_PROGRAMS=OFF \
            -DZSTD_BUILD_SHARED=OFF \
            -DZSTD_BUILD_STATIC=ON \
            -DZSTD_BUILD_TESTS=OFF \
            -DZSTD_MULTITHREAD_SUPPORT=OFF \
            -DZSTD_LEGACY_SUPPORT=OFF
        cmake --build "$WORK_DIR/build-zstd-$flavour" -j "$JOBS"
        cmake --install "$WORK_DIR/build-zstd-$flavour"
    fi

    if [ "$FORCE" = "no" ] && [ -f "$fprefix/lib/liblz4.a" ]; then
        echo "== lz4 $LZ4_VERSION ($flavour) already installed, skipping =="
    else
        echo "== lz4 $LZ4_VERSION ($flavour) =="
        fetch "https://github.com/lz4/lz4/releases/download/v$LZ4_VERSION/lz4-$LZ4_VERSION.tar.gz" \
              "lz4-$LZ4_VERSION.tar.gz" "$LZ4_SHA256" || exit 1
        rm -rf "$WORK_DIR/lz4-$LZ4_VERSION" "$WORK_DIR/build-lz4-$flavour"
        tar -xzf "$DL_DIR/lz4-$LZ4_VERSION.tar.gz" -C "$WORK_DIR"
        CFLAGS="-Oz $pt" emcmake cmake -S "$WORK_DIR/lz4-$LZ4_VERSION/build/cmake" \
            -B "$WORK_DIR/build-lz4-$flavour" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="$fprefix" \
            -DLZ4_BUILD_CLI=OFF \
            -DLZ4_BUILD_LEGACY_LZ4C=OFF \
            -DBUILD_SHARED_LIBS=OFF \
            -DBUILD_STATIC_LIBS=ON
        cmake --build "$WORK_DIR/build-lz4-$flavour" -j "$JOBS"
        cmake --install "$WORK_DIR/build-lz4-$flavour"
    fi

    # ------------------------------------------------------------------
    # KaHIP (the serial kaffpa partitioner, full artifact only)
    #
    # Only the static `kahip_static` interface is linked (PARHIP off, no MPI
    # -- the same link surface as the native KaHIP integration, which links
    # no ParHIP either). OpenMP is off under Emscripten even for the mt
    # flavour: the serial kaffpa path is what ships (deterministic under a
    # fixed seed), and Emscripten's OpenMP needs -fopenmp on every unit,
    # which the $pt-only mt flags do not provide.
    # ------------------------------------------------------------------
    if [ "$FORCE" = "no" ] && [ -f "$fprefix/lib/libkahip_static.a" ]; then
        echo "== KaHIP $KAHIP_VERSION ($flavour) already installed, skipping =="
    else
        echo "== KaHIP $KAHIP_VERSION ($flavour) =="
        rm -rf "$WORK_DIR/KaHIP-$KAHIP_VERSION" "$WORK_DIR/build-kahip-$flavour"
        unpack_tag "https://github.com/KaHIP/KaHIP/archive/refs/tags/v$KAHIP_VERSION.tar.gz" \
              "kahip-$KAHIP_VERSION.tar.gz" "$KAHIP_SHA256" \
              "https://github.com/KaHIP/KaHIP.git" "v$KAHIP_VERSION" "$KAHIP_COMMIT" \
              "KaHIP-$KAHIP_VERSION"
        apply_patches "kahip-$KAHIP_VERSION" "$WORK_DIR/KaHIP-$KAHIP_VERSION"
        CFLAGS="-Oz $pt" CXXFLAGS="-Oz $pt" emcmake cmake -S "$WORK_DIR/KaHIP-$KAHIP_VERSION" \
            -B "$WORK_DIR/build-kahip-$flavour" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="$fprefix" \
            -DPARHIP=OFF -DNOMPI=ON -DNONATIVEOPTIMIZATIONS=ON \
            -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE
        cmake --build "$WORK_DIR/build-kahip-$flavour" -j "$JOBS"
        cmake --install "$WORK_DIR/build-kahip-$flavour"
    fi

    # ------------------------------------------------------------------
    # Netgen (nglib and ngcore, the tetrahedralize backend, full artifact only)
    #
    # Static, -fwasm-exceptions like the core (nglib throws ngcore::Exception
    # through Ng_GenerateVolumeMesh). Two adjustments, both for Emscripten
    # 6.0.3: the patch drops the JS exception catching and the web-only link
    # settings ngcore imposes on its consumers, and -DEMSCRIPTEN restores the
    # macro its sources still test (Emscripten predefines only
    # __EMSCRIPTEN__ now), which keeps the thread-priority calls the seq
    # flavour has no symbols for out of the task manager. The version is
    # passed in: the sources are a git archive, which has no history for
    # Netgen's `git describe`.
    # ------------------------------------------------------------------
    if [ "$FORCE" = "no" ] && [ -f "$fprefix/lib/libnglib.a" ]; then
        echo "== Netgen $NETGEN_VERSION ($flavour) already installed, skipping =="
    else
        echo "== Netgen $NETGEN_VERSION ($flavour) =="
        fetch_git "https://github.com/NGSolve/netgen.git" "v$NETGEN_VERSION" "$NETGEN_COMMIT" \
            "netgen-$NETGEN_VERSION.git"
        rm -rf "$WORK_DIR/netgen-$NETGEN_VERSION" "$WORK_DIR/build-netgen-$flavour"
        mkdir -p "$WORK_DIR/netgen-$NETGEN_VERSION"
        (cd "$DL_DIR/netgen-$NETGEN_VERSION.git" && git archive --format=tar HEAD) |
            tar -xf - -C "$WORK_DIR/netgen-$NETGEN_VERSION"
        apply_patches "netgen-$NETGEN_VERSION" "$WORK_DIR/netgen-$NETGEN_VERSION"
        CFLAGS="-DEMSCRIPTEN -sUSE_ZLIB=1 $pt" \
        CXXFLAGS="-DEMSCRIPTEN -fwasm-exceptions -sUSE_ZLIB=1 $pt" \
        emcmake cmake -S "$WORK_DIR/netgen-$NETGEN_VERSION" \
            -B "$WORK_DIR/build-netgen-$flavour" \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="$fprefix" \
            -DNETGEN_VERSION_GIT="v$NETGEN_VERSION-0-g$(echo "$NETGEN_COMMIT" | cut -c1-7)" \
            -DUSE_SUPERBUILD=OFF -DUSE_GUI=OFF -DUSE_PYTHON=OFF -DUSE_OCC=OFF \
            -DUSE_MPI=OFF -DUSE_NATIVE_ARCH=OFF -DUSE_JPEG=OFF -DUSE_MPEG=OFF \
            -DUSE_CGNS=OFF -DBUILD_STUB_FILES=OFF -DENABLE_UNIT_TESTS=OFF \
            -DNGLIB_LIBRARY_TYPE=STATIC -DNGCORE_LIBRARY_TYPE=STATIC \
            -DNGGUI_LIBRARY_TYPE=STATIC \
            -DZLIB_INCLUDE_DIRS="$SYSROOT/include"
        cmake --build "$WORK_DIR/build-netgen-$flavour" -j "$JOBS"
        cmake --install "$WORK_DIR/build-netgen-$flavour"
        # Its tutorials and manual are installed too; the link needs none.
        rm -rf "$WORK_DIR/build-netgen-$flavour" "$fprefix/share/netgen"
    fi

    [ "$WITH_OCCT" = "yes" ] || return 0
    if [ "$FORCE" = "no" ] && [ -f "$fprefix/lib/libTKernel.a" ]; then
        echo "== OpenCASCADE $OCCT_VERSION ($flavour) already installed, skipping =="
        return 0
    fi
    echo "== OpenCASCADE $OCCT_VERSION ($flavour) =="
    rm -rf "$WORK_DIR/build-occt-$flavour"
    occt_source "$DL_DIR" "$WORK_DIR/occt-$OCCT_VERSION"
    # OCCT defines OCC_CONVERT_SIGNALS on every non-MSVC build: OCC_CATCH_SIGNALS
    # then setjmp()s so a signal handler can longjmp a SIGSEGV or SIGFPE into a
    # C++ exception. WebAssembly raises no signals, and the setjmp/longjmp it
    # lowers to (Wasm SjLj, a second exception tag) is code Binaryen 6.0.3's
    # wasm-opt fails to parse at link ("popping from empty stack"). The
    # flags' -U follows CMake's -D on every compile line, so OCC_CATCH_SIGNALS
    # is empty and Standard_Failure is a plain C++ exception throughout.
    CFLAGS="$pt" CXXFLAGS="$pt" OCCT_CMAKE="emcmake cmake" occt_configure \
        "$WORK_DIR/occt-$OCCT_VERSION" "$WORK_DIR/build-occt-$flavour" "$fprefix" Static \
        -DCMAKE_CXX_FLAGS="-fwasm-exceptions -UOCC_CONVERT_SIGNALS $pt" \
        -DCMAKE_C_FLAGS="-UOCC_CONVERT_SIGNALS $pt"
    cmake --build "$WORK_DIR/build-occt-$flavour" -j "$JOBS"
    cmake --install "$WORK_DIR/build-occt-$flavour"
    # Nor do its consumers get it: the installed config would add the define
    # to every directory that finds OpenCASCADE.
    sed -i.bak '/OCC_CONVERT_SIGNALS/d' \
        "$fprefix/lib/cmake/opencascade/OpenCASCADECompileDefinitionsAndFlags-release.cmake"
    rm -f "$fprefix/lib/cmake/opencascade/OpenCASCADECompileDefinitionsAndFlags-release.cmake.bak"
    rm -rf "$WORK_DIR/build-occt-$flavour"
}

if [ "$WITH_FULL" = "yes" ]; then
    # The lean section's CFLAGS (-Oz, the zlib port) are not the full
    # libraries' flags; each build above sets its own.
    unset CFLAGS CXXFLAGS
    for flavour in $FULL_FLAVOURS; do
        build_full_flavour "$flavour"
    done
fi

# The unpacked sources and CMake build trees are ~1 GB and re-derivable from
# the (checksummed) tarballs in _dl, so they are not kept -- what matters, and
# what CI caches, is the install prefix.
rm -rf "$WORK_DIR"

echo
echo "== done =="
if [ "$WITH_LEAN" = "yes" ]; then
    ls -l "$PREFIX/lib/"*.a 2>/dev/null || true
    echo "prefix: $PREFIX"
fi
if [ "$WITH_FULL" = "yes" ]; then
    for flavour in $FULL_FLAVOURS; do
        echo "full prefix ($flavour): $(full_prefix "$flavour")"
    done
fi
echo "next:   $SCRIPT_DIR/configure-wasm.sh --build"
