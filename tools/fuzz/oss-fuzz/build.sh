#!/bin/bash -eu
# OSS-Fuzz build script (doc/fuzzing.md): one fuzz target per native
# reader, as copies of the single harness named meshioplusplus_fuzz_read_<fmt>
# (the harness takes its format from its own name). $CC/$CXX/$CFLAGS and
# $LIB_FUZZING_ENGINE come from the OSS-Fuzz environment. No local sanitizer
# defaults override the sanitizer selected by OSS-Fuzz.
BUILD="$WORK/meshioplusplus-build"
cmake -S "$SRC/meshioplusplus" -B "$BUILD" -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF -DMESHIOPLUSPLUS_BUILD_TESTS=OFF \
    -DMESHIOPLUSPLUS_BUILD_FUZZERS=ON \
    -DMESHIOPLUSPLUS_FUZZING_ENGINE="$LIB_FUZZING_ENGINE" \
    -DMESHIOPLUSPLUS_WITH_HDF5=OFF -DMESHIOPLUSPLUS_WITH_NETCDF=OFF \
    -DMESHIOPLUSPLUS_WITH_GIDPOST=OFF -DMESHIOPLUSPLUS_WITH_BZIP2=OFF \
    -DMESHIOPLUSPLUS_ZLIB_STATIC=ON -DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ
cmake --build "$BUILD" --parallel "${JOBS:-4}" \
    --target meshioplusplus_fuzz_read meshioplusplus_fuzz_replay meshioplusplus_fuzz_seeds
while IFS= read -r fmt; do
    grep -qx "$fmt" "$SRC/meshioplusplus/tools/fuzz/not_fuzzed.txt" && continue
    cp "$BUILD/meshioplusplus_fuzz_read" "$OUT/meshioplusplus_fuzz_read_$fmt"
    [ -f "$SRC/meshioplusplus/tools/fuzz/dicts/$fmt.dict" ] &&
        cp "$SRC/meshioplusplus/tools/fuzz/dicts/$fmt.dict" "$OUT/meshioplusplus_fuzz_read_$fmt.dict"
    # A format without a dictionary must not leave the loop with status 1.
    true
done < <("$BUILD/meshioplusplus_fuzz_replay" -list-formats)
"$BUILD/meshioplusplus_fuzz_seeds" "$WORK/generated-seeds"
python3 "$SRC/meshioplusplus/tools/fuzz/oss-fuzz/package_seeds.py" \
    "$SRC/meshioplusplus/tests/fuzz/regressions" "$OUT" --generated "$WORK/generated-seeds"
