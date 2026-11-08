#!/bin/bash
# Native, instrumented HDF5/netCDF/zlib prefix for the library-reader campaign.
# Unlike the core's CMake, this opt-in script downloads pinned dependencies.
# Usage: CC=clang CXX=clang++ tools/fuzz/build_library_deps.sh PREFIX [WORK_DIR]
set -euo pipefail
PREFIX=$(realpath -m "${1:?usage: $0 PREFIX [WORK_DIR]}")
WORK=$(realpath -m "${2:-$PREFIX-work}")
export CC=${CC:-clang} CXX=${CXX:-clang++}
export CFLAGS=${CFLAGS:--O1 -g -fno-omit-frame-pointer -fsanitize=fuzzer-no-link,address,undefined -fno-sanitize-recover=undefined}
export CXXFLAGS=${CXXFLAGS:-$CFLAGS}
JOBS=${JOBS:-4}
mkdir -p "$PREFIX" "$WORK"
patch_hash=$(sha256sum "$(dirname "$0")/patches/netcdf-4.9.3-utf8proc-null-guard.patch" | cut -d' ' -f1)
manifest=$(printf '%s\n' "recipe=4 zlib=1.3.1 hdf5=1.14.6 netcdf=4.9.3 hdf5-no-function-check netcdf-null-guard=$patch_hash" \
    "$($CC --version)" "$($CXX --version)" "$CFLAGS" "$CXXFLAGS")
if [ -f "$PREFIX/instrumentation.txt" ]; then
    [ "$(cat "$PREFIX/instrumentation.txt")" = "$manifest" ] || {
        echo "Compiler/flags changed: use a fresh dependency prefix." >&2; exit 1;
    }
    exit 0
fi
# fetch <archive> <sha256> <url>...: the first URL that serves the archive; the
# checksum is checked whichever served it.
fetch() {
    local archive=$1 hash=$2 url
    shift 2
    if [ ! -f "$WORK/$archive" ]; then
        for url in "$@"; do
            curl --fail --location --retry 3 "$url" -o "$WORK/$archive" && break
            rm -f "$WORK/$archive"
        done
        [ -f "$WORK/$archive" ] || { echo "error: no URL served $archive" >&2; exit 1; }
    fi
    printf '%s  %s\n' "$hash" "$WORK/$archive" | sha256sum --check
    tar -xf "$WORK/$archive" -C "$WORK"
}
# fetch_tag <archive> <sha256> <url> <git-url> <tag> <commit> <dir>: the
# checksummed archive, or, when no URL serves it, a shallow clone of <tag>
# refused unless it is <commit> (the pins of build/deps-pins.sh). GitHub's
# release storage has gone down on its own while clones kept working (504s
# for hours, November 2026), and some proxies refuse codeload's tag archives.
fetch_tag() {
    if ! (fetch "$1" "$2" "$3"); then
        echo "$1 unavailable, cloning $4 $5" >&2
        rm -rf "${WORK:?}/$7"
        git -c advice.detachedHead=false clone --quiet --depth 1 --branch "$5" "$4" "$WORK/$7"
        [ "$(git -C "$WORK/$7" rev-parse HEAD)" = "$6" ] || {
            echo "error: $4 $5 is not the pinned commit $6" >&2; exit 1;
        }
    fi
}
# zlib has zlib.net's archive of the same tarball.
fetch zlib-1.3.1.tar.gz 9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23 \
    https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz \
    https://zlib.net/fossils/zlib-1.3.1.tar.gz
fetch_tag hdf5-1.14.6.tar.gz e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b \
    https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz \
    https://github.com/HDFGroup/hdf5.git hdf5_1.14.6 7bf340440909d468dbb3cf41f0ea0d87f5050cea \
    hdf5-1.14.6
fetch_tag netcdf-c-4.9.3.tar.gz 990f46d49525d6ab5dc4249f8684c6deeaf54de6fec63a187e9fb382cc0ffdff \
    https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.3.tar.gz \
    https://github.com/Unidata/netcdf-c.git v4.9.3 3a6b1debf1557f07b606ce3653e44f0d711203be \
    netcdf-c-4.9.3
# Backported upstream maintenance fix for the bundled Unicode helper's
# measuring pass (see tools/fuzz/patches/). Applied after every extract so
# a reused work directory cannot keep an unpatched tree.
patch -d "$WORK/netcdf-c-4.9.3" -N -p1 --no-backup-if-mismatch \
    < "$(dirname "$0")/patches/netcdf-4.9.3-utf8proc-null-guard.patch"
common=(-G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="$PREFIX"
    -DCMAKE_PREFIX_PATH="$PREFIX" -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=ON
    -DCMAKE_INSTALL_RPATH="$PREFIX/lib")
cmake -S "$WORK/zlib-1.3.1" -B "$WORK/build-zlib" "${common[@]}" -DZLIB_BUILD_EXAMPLES=OFF
cmake --build "$WORK/build-zlib" --parallel "$JOBS"
cmake --install "$WORK/build-zlib"
cmake -S "$WORK/hdf5-1.14.6" -B "$WORK/build-hdf5" "${common[@]}" \
    -DCMAKE_C_FLAGS="$CFLAGS -fno-sanitize=function" \
    -DBUILD_TESTING=OFF -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_UTILS=OFF \
    -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_FORTRAN=OFF \
    -DHDF5_BUILD_JAVA=OFF -DHDF5_BUILD_HL_LIB=ON -DHDF5_ENABLE_PARALLEL=OFF \
    -DHDF5_ENABLE_PLUGIN_SUPPORT=OFF -DHDF5_ENABLE_ROS3_VFD=OFF \
    -DHDF5_ENABLE_SZIP_SUPPORT=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT=ON \
    -DZLIB_INCLUDE_DIR="$PREFIX/include" -DZLIB_LIBRARY="$PREFIX/lib/libz.so"
cmake --build "$WORK/build-hdf5" --parallel "$JOBS"
cmake --install "$WORK/build-hdf5"
cmake -S "$WORK/netcdf-c-4.9.3" -B "$WORK/build-netcdf" "${common[@]}" \
    -DNETCDF_ENABLE_HDF5=ON -DNETCDF_ENABLE_DAP=OFF -DNETCDF_ENABLE_DAP4=OFF \
    -DNETCDF_ENABLE_BYTERANGE=OFF -DNETCDF_ENABLE_NCZARR=OFF -DNETCDF_ENABLE_S3=OFF \
    -DNETCDF_ENABLE_LIBXML2=OFF -DNETCDF_ENABLE_PLUGINS=OFF -DNETCDF_ENABLE_TESTS=OFF \
    -DNETCDF_ENABLE_FILTER_SZIP=OFF -DNETCDF_ENABLE_FILTER_BZ2=OFF \
    -DNETCDF_ENABLE_FILTER_BLOSC=OFF -DNETCDF_ENABLE_FILTER_ZSTD=OFF \
    -DNETCDF_ENABLE_EXAMPLES=OFF -DNETCDF_BUILD_UTILITIES=OFF -DBUILD_TESTING=OFF
cmake --build "$WORK/build-netcdf" --parallel "$JOBS"
cmake --install "$WORK/build-netcdf"
printf '%s\n' "$manifest" > "$PREFIX/instrumentation.txt"
