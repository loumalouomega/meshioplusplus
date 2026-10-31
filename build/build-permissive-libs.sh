#!/bin/sh
# build-permissive-libs.sh -- build the permissively licensed libraries the PyPI
# base wheel links natively (roadmap 4.7.2) as static, position-independent
# archives in one prefix: zlib is vendored by the main CMake, so this builds
# bzip2, zstd, lz4, HDF5, netCDF, cgnslib and KaHIP. Eigen is a header-only
# submodule and needs nothing here.
#
#   build/build-permissive-libs.sh PREFIX
#
# Then configure meshio++ with
#   -DMESHIOPLUSPLUS_PREFER_STATIC_DEPS=ON -DCMAKE_PREFIX_PATH=PREFIX
#   -DMESHIOPLUSPLUS_WITH_HDF5=ON -DMESHIOPLUSPLUS_WITH_NETCDF=ON
#   -DMESHIOPLUSPLUS_WITH_CGNSLIB=ON -DMESHIOPLUSPLUS_WITH_ZSTD=ON
#   -DMESHIOPLUSPLUS_WITH_LZ4=ON -DMESHIOPLUSPLUS_WITH_BZIP2=ON
#   -DMESHIOPLUSPLUS_WITH_KAHIP=ON -DKAHIP_ROOT=PREFIX
#
# Versions and checksums are build/deps-pins.sh's, shared with the WebAssembly
# build (build-wasm-deps.sh). Every library is built with -fvisibility=hidden
# and position-independent code: the wheel's shared core links the archives and
# keeps their symbols local (MESHIOPLUSPLUS_HIDDEN_ARCHIVES), because h5py and
# netCDF4 bring their own HDF5 into the same process.
#
# Like build-occt.sh, a shell script and not CMake: meshio++'s CMake never
# downloads anything, it only finds an already-built dependency. The prefix is
# not relocatable (netCDF's installed config records the HDF5 path).

set -eu

case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*)
        # The flags below are GCC/Clang's, and bzip2 builds with make; the Windows
        # wheels keep h5py and netCDF4 for these formats until an MSVC recipe exists.
        echo "build-permissive-libs: Windows (MSVC) is not supported yet" >&2
        exit 2 ;;
esac

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
. "$SCRIPT_DIR/deps-pins.sh"
BZIP2_VERSION="1.0.8"
BZIP2_SHA256="ab5a03176ee106d3f0fa90e381da478ddae405918153cca248e682cd0c4a2269"

PREFIX=""
UNIVERSAL="no"
FORCE="no"
JOBS=""
WORK=""
ONLY=""

usage() {
    cat <<EOF
Usage: build/build-permissive-libs.sh PREFIX [options]
  --universal      macOS: build x86_64 and arm64 (deployment target 13.3)
  --only LIB       build one of: bzip2 zstd lz4 hdf5 netcdf cgns kahip
  --work DIR       download and build directory (default: PREFIX/_work)
  --jobs N         parallel build jobs (default: all cores)
  --force          rebuild even when PREFIX already holds the library
Builds static, hidden-visibility archives: bzip2 $BZIP2_VERSION, zstd $ZSTD_VERSION,
lz4 $LZ4_VERSION, HDF5 $HDF5_VERSION, netCDF $NETCDF_VERSION, cgnslib $CGNS_VERSION, KaHIP $KAHIP_VERSION.
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --universal) UNIVERSAL="yes"; shift ;;
        --only) ONLY="$2"; shift 2 ;;
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
DL="$WORK/downloads"
mkdir -p "$WORK" "$DL"
[ -n "$JOBS" ] || JOBS=$( (command -v nproc >/dev/null 2>&1 && nproc) ||
                          (sysctl -n hw.ncpu 2>/dev/null) || echo 2 )

wanted() { [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }
have() { [ "$FORCE" = "no" ] && [ -f "$PREFIX/lib/$1" ]; }

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

# fetch URL FILE SHA256: a verified download into the cache.
fetch() {
    if [ ! -f "$DL/$2" ] || [ "$(sha256_of "$DL/$2")" != "$3" ]; then
        echo "  downloading $2"
        curl -fsSL --retry 3 -o "$DL/$2" "$1"
    fi
    got=$(sha256_of "$DL/$2")
    if [ "$got" != "$3" ]; then
        echo "error: $2 is $got, expected $3" >&2
        rm -f "$DL/$2"
        exit 1
    fi
}

# unpack_tag ARCHIVE_URL FILE SHA256 GIT_URL TAG COMMIT DIR: the release
# archive, or -- where GitHub's codeload is refused -- a shallow clone of the
# tag verified against the commit it must peel to.
unpack_tag() {
    rm -rf "$WORK/$7"
    if (fetch "$1" "$2" "$3") 2>/dev/null; then
        mkdir -p "$WORK/$7"
        tar -xzf "$DL/$2" -C "$WORK/$7" --strip-components=1
    else
        git -c advice.detachedHead=false clone --quiet --depth 1 --branch "$5" "$4" "$WORK/$7"
        got=$(git -C "$WORK/$7" rev-parse HEAD)
        if [ "$got" != "$6" ]; then
            echo "error: $4 $5 is $got, expected $6" >&2
            exit 1
        fi
    fi
}

# Every archive is position-independent with hidden symbols; the one macOS
# flavour builds both architectures.
CFLAGS_ALL="-fPIC -fvisibility=hidden -O2"
CXXFLAGS_ALL="-fPIC -fvisibility=hidden -O2"
CMAKE_COMMON="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_PREFIX_PATH=$PREFIX
 -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_C_VISIBILITY_PRESET=hidden
 -DCMAKE_CXX_VISIBILITY_PRESET=hidden -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DBUILD_TESTING=OFF"
case "$(uname -s)" in
    Darwin)
        CMAKE_COMMON="$CMAKE_COMMON -DCMAKE_OSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET:-13.3}"
        [ "$UNIVERSAL" = "no" ] || CMAKE_COMMON="$CMAKE_COMMON -DCMAKE_OSX_ARCHITECTURES=x86_64;arm64"
        ;;
esac
# shellcheck disable=SC2086  # CMAKE_COMMON is a word list on purpose
cmake_build() {  # cmake_build NAME SRC [extra args]: configure, build, install
    name="$1"; src="$2"; shift 2
    rm -rf "$WORK/build-$name"
    CFLAGS="$CFLAGS_ALL" CXXFLAGS="$CXXFLAGS_ALL" cmake -S "$src" -B "$WORK/build-$name" $CMAKE_COMMON "$@"
    cmake --build "$WORK/build-$name" --config Release -j "$JOBS"
    cmake --install "$WORK/build-$name" --config Release
    rm -rf "$WORK/build-$name"
}

# -- bzip2 (no CMake in 1.0.8: its Makefile builds the archive) ---------------
if wanted bzip2; then
    if have libbz2.a; then echo "== bzip2 already installed =="; else
        echo "== bzip2 $BZIP2_VERSION =="
        fetch "https://sourceware.org/pub/bzip2/bzip2-$BZIP2_VERSION.tar.gz" \
              "bzip2-$BZIP2_VERSION.tar.gz" "$BZIP2_SHA256"
        rm -rf "$WORK/bzip2"; mkdir -p "$WORK/bzip2"
        tar -xzf "$DL/bzip2-$BZIP2_VERSION.tar.gz" -C "$WORK/bzip2" --strip-components=1
        (cd "$WORK/bzip2" &&
         make -j "$JOBS" libbz2.a CFLAGS="$CFLAGS_ALL -D_FILE_OFFSET_BITS=64" &&
         mkdir -p "$PREFIX/lib" "$PREFIX/include" &&
         cp libbz2.a "$PREFIX/lib/" && cp bzlib.h "$PREFIX/include/")
    fi
fi

# -- zstd ----------------------------------------------------------------------
if wanted zstd; then
    if have libzstd.a; then echo "== zstd already installed =="; else
        echo "== zstd $ZSTD_VERSION =="
        fetch "https://github.com/facebook/zstd/releases/download/v$ZSTD_VERSION/zstd-$ZSTD_VERSION.tar.gz" \
              "zstd-$ZSTD_VERSION.tar.gz" "$ZSTD_SHA256"
        rm -rf "$WORK/zstd"; mkdir -p "$WORK/zstd"
        tar -xzf "$DL/zstd-$ZSTD_VERSION.tar.gz" -C "$WORK/zstd" --strip-components=1
        cmake_build zstd "$WORK/zstd/build/cmake" \
            -DZSTD_BUILD_PROGRAMS=OFF -DZSTD_BUILD_SHARED=OFF -DZSTD_BUILD_STATIC=ON \
            -DZSTD_BUILD_TESTS=OFF -DZSTD_MULTITHREAD_SUPPORT=OFF -DZSTD_LEGACY_SUPPORT=OFF
    fi
fi

# -- lz4 -----------------------------------------------------------------------
if wanted lz4; then
    if have liblz4.a; then echo "== lz4 already installed =="; else
        echo "== lz4 $LZ4_VERSION =="
        fetch "https://github.com/lz4/lz4/releases/download/v$LZ4_VERSION/lz4-$LZ4_VERSION.tar.gz" \
              "lz4-$LZ4_VERSION.tar.gz" "$LZ4_SHA256"
        rm -rf "$WORK/lz4"; mkdir -p "$WORK/lz4"
        tar -xzf "$DL/lz4-$LZ4_VERSION.tar.gz" -C "$WORK/lz4" --strip-components=1
        cmake_build lz4 "$WORK/lz4/build/cmake" \
            -DLZ4_BUILD_CLI=OFF -DLZ4_BUILD_LEGACY_LZ4C=OFF \
            -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON
    fi
fi

# -- HDF5 ----------------------------------------------------------------------
# Serial, no tools, with its own zlib found on the system (the wheel's static
# zlib is vendored by the main CMake and linked there). h5py ships its own copy.
if wanted hdf5; then
    if have libhdf5.a; then echo "== HDF5 already installed =="; else
        echo "== HDF5 $HDF5_VERSION =="
        fetch "https://github.com/HDFGroup/hdf5/releases/download/hdf5_$HDF5_VERSION/hdf5-$HDF5_VERSION.tar.gz" \
              "hdf5-$HDF5_VERSION.tar.gz" "$HDF5_SHA256"
        # The release archive's layout is not the usual single top directory
        # (`--strip-components=1` left no CMakeLists.txt), so find the source
        # root instead of assuming its depth.
        rm -rf "$WORK/hdf5" "$WORK/hdf5-unpack"; mkdir -p "$WORK/hdf5-unpack"
        tar -xzf "$DL/hdf5-$HDF5_VERSION.tar.gz" -C "$WORK/hdf5-unpack"
        root=$(find "$WORK/hdf5-unpack" -maxdepth 3 -name CMakeLists.txt -print | awk '{ print length, $0 }' | sort -n | head -n 1 | cut -d' ' -f2-)
        [ -n "$root" ] || { echo "error: no CMakeLists.txt in the HDF5 archive" >&2; exit 1; }
        mv "$(dirname "$root")" "$WORK/hdf5"
        rm -rf "$WORK/hdf5-unpack"
        cmake_build hdf5 "$WORK/hdf5" \
            -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DONLY_SHARED_LIBS=OFF \
            -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_UTILS=OFF -DHDF5_BUILD_EXAMPLES=OFF \
            -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_FORTRAN=OFF -DHDF5_BUILD_JAVA=OFF \
            -DHDF5_BUILD_HL_LIB=ON -DHDF5_ENABLE_DEPRECATED_SYMBOLS=ON \
            -DHDF5_ENABLE_THREADSAFE=OFF -DHDF5_ENABLE_PARALLEL=OFF \
            -DHDF5_ENABLE_PLUGIN_SUPPORT=OFF -DHDF5_ENABLE_ROS3_VFD=OFF \
            -DHDF5_ENABLE_DIRECT_VFD=OFF -DHDF5_ENABLE_MIRROR_VFD=OFF \
            -DHDF5_ENABLE_SUBFILING_VFD=OFF -DHDF5_ENABLE_SZIP_SUPPORT=OFF \
            -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DHDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16=OFF
    fi
fi

# -- netCDF-C ------------------------------------------------------------------
# Classic netCDF-3 plus netCDF-4 over the static HDF5 above: what exodus.cpp
# needs. Every option under both spellings, as in build-wasm-deps.sh (4.9.3
# renamed them, and an unknown -D is only a warning).
if wanted netcdf; then
    if have libnetcdf.a; then echo "== netCDF already installed =="; else
        echo "== netCDF $NETCDF_VERSION =="
        unpack_tag "https://github.com/Unidata/netcdf-c/archive/refs/tags/v$NETCDF_VERSION.tar.gz" \
                   "netcdf-c-$NETCDF_VERSION.tar.gz" "$NETCDF_SHA256" \
                   "https://github.com/Unidata/netcdf-c.git" "v$NETCDF_VERSION" "$NETCDF_COMMIT" netcdf-c
        cmake_build netcdf "$WORK/netcdf-c" \
            -DBUILD_SHARED_LIBS=OFF -DHDF5_USE_STATIC_LIBRARIES=ON \
            -DENABLE_HDF5=ON -DNETCDF_ENABLE_HDF5=ON \
            -DENABLE_DAP=OFF -DNETCDF_ENABLE_DAP=OFF \
            -DENABLE_DAP4=OFF -DNETCDF_ENABLE_DAP4=OFF \
            -DENABLE_BYTERANGE=OFF -DNETCDF_ENABLE_BYTERANGE=OFF \
            -DENABLE_NCZARR=OFF -DNETCDF_ENABLE_NCZARR=OFF \
            -DENABLE_S3=OFF -DNETCDF_ENABLE_S3=OFF \
            -DENABLE_LIBXML2=OFF -DNETCDF_ENABLE_LIBXML2=OFF \
            -DENABLE_PLUGINS=OFF -DNETCDF_ENABLE_PLUGINS=OFF \
            -DENABLE_FILTER_TESTING=OFF -DNETCDF_ENABLE_FILTER_TESTING=OFF \
            -DENABLE_TESTS=OFF -DNETCDF_ENABLE_TESTS=OFF \
            -DENABLE_EXAMPLES=OFF -DNETCDF_ENABLE_EXAMPLES=OFF \
            -DBUILD_UTILITIES=OFF -DNETCDF_BUILD_UTILITIES=OFF
    fi
fi

# -- cgnslib -------------------------------------------------------------------
# The top-level source directory, not src/: the root CMakeLists runs the
# CHECK_TYPE_SIZE calls that src/CMakeLists.txt reads to pick cgsize_t.
if wanted cgns; then
    if have libcgns.a; then echo "== cgnslib already installed =="; else
        echo "== cgnslib $CGNS_VERSION =="
        unpack_tag "https://github.com/CGNS/CGNS/archive/refs/tags/v$CGNS_VERSION.tar.gz" \
                   "cgns-$CGNS_VERSION.tar.gz" "$CGNS_SHA256" \
                   "https://github.com/CGNS/CGNS.git" "v$CGNS_VERSION" "$CGNS_COMMIT" CGNS
        cmake_build cgns "$WORK/CGNS" \
            -DBUILD_SHARED_LIBS=OFF -DCGNS_BUILD_SHARED=OFF -DCGNS_USE_SHARED=OFF \
            -DCGNS_ENABLE_HDF5=ON -DCGNS_ENABLE_FORTRAN=OFF -DCGNS_ENABLE_PARALLEL=OFF \
            -DCGNS_ENABLE_TESTS=OFF -DCGNS_BUILD_TESTING=OFF -DCGNS_BUILD_CGNSTOOLS=OFF \
            -DCGNS_ENABLE_64BIT=ON -DHDF5_NEED_ZLIB=ON \
            -DHDF5_INCLUDE_PATH="$PREFIX/include" -DHDF5_LIBRARY="$PREFIX/lib/libhdf5.a"
    fi
fi

# -- KaHIP (the serial kaffpa interface: no ParHIP, no MPI, no OpenMP) ----------
if wanted kahip; then
    if have libkahip_static.a || have libkahip.a; then echo "== KaHIP already installed =="; else
        echo "== KaHIP $KAHIP_VERSION =="
        unpack_tag "https://github.com/KaHIP/KaHIP/archive/refs/tags/v$KAHIP_VERSION.tar.gz" \
                   "kahip-$KAHIP_VERSION.tar.gz" "$KAHIP_SHA256" \
                   "https://github.com/KaHIP/KaHIP.git" "v$KAHIP_VERSION" "$KAHIP_COMMIT" KaHIP
        cmake_build kahip "$WORK/KaHIP" \
            -DPARHIP=OFF -DNOMPI=ON -DNONATIVEOPTIMIZATIONS=ON \
            -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE
    fi
fi

echo "== done: -DMESHIOPLUSPLUS_PREFER_STATIC_DEPS=ON -DCMAKE_PREFIX_PATH=$PREFIX -DKAHIP_ROOT=$PREFIX =="
