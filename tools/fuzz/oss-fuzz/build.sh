#!/bin/bash -eu
# Draft OSS-Fuzz build script (doc/fuzzing.md): one fuzz target per native
# reader, as copies of the single harness named meshioplusplus_fuzz_read_<fmt>
# (the harness takes its format from its own name). $CC/$CXX/$CFLAGS and
# $LIB_FUZZING_ENGINE come from the OSS-Fuzz environment; the HDF5/netCDF
# formats are left out (tools/fuzz/not_fuzzed.txt).
cmake -S "$SRC/meshioplusplus" -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF -DMESHIOPLUSPLUS_BUILD_TESTS=ON \
    -DMESHIOPLUSPLUS_WITH_HDF5=OFF -DMESHIOPLUSPLUS_WITH_NETCDF=OFF \
    -DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ
cmake --build build --target meshioplusplus_core_obj meshioplusplus_fuzz_replay
$CXX $CXXFLAGS -std=c++20 -I"$SRC/meshioplusplus/src/cpp/include" \
    "$SRC/meshioplusplus/tests/fuzz/fuzz_read.cpp" \
    build/CMakeFiles/meshioplusplus_core_obj.dir/src/cpp/src/*.o \
    build/CMakeFiles/meshioplusplus_core_obj.dir/src/cpp/src/*/*.o \
    $LIB_FUZZING_ENGINE -o "$OUT/meshioplusplus_fuzz_read"
skip=$(grep -v '^#' "$SRC/meshioplusplus/tools/fuzz/not_fuzzed.txt" | grep -v '^$')
for fmt in $(build/meshioplusplus_fuzz_replay -list-formats); do
    echo "$skip" | grep -qx "$fmt" && continue
    cp "$OUT/meshioplusplus_fuzz_read" "$OUT/meshioplusplus_fuzz_read_$fmt"
    [ -f "$SRC/meshioplusplus/tools/fuzz/dicts/$fmt.dict" ] &&
        cp "$SRC/meshioplusplus/tools/fuzz/dicts/$fmt.dict" "$OUT/meshioplusplus_fuzz_read_$fmt.dict"
done
rm "$OUT/meshioplusplus_fuzz_read"
