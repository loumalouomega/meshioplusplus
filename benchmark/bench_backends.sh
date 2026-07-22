#!/bin/sh
# bench_backends.sh — build and run the C++ mesh-backend benchmark
# (src/cpp/benchmark/bench_backends.cpp) for every mesh backend and collate the
# results into benchmark/results_backends.csv.
#
# The mesh backend (MESHIO / NATIVE / KRATOS) is an exclusive compile-time
# choice, so one build tree per backend is configured under
# build/bench-<backend> (Python extension off, benchmark target on, parallel
# backend fixed to OPENMP so only the mesh backend varies across runs).
#
#   ./benchmark/bench_backends.sh          # default grid size (n=35, 257k tets)
#   ./benchmark/bench_backends.sh 50       # bigger grid
#
# Output columns: backend,op,format,cells,median_s,runs

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_DIR=$(dirname -- "$SCRIPT_DIR")
N="${1:-35}"
OUT="$SCRIPT_DIR/results_backends.csv"
PARALLEL_BACKEND="${MESHIOPLUSPLUS_BENCH_PARALLEL:-OPENMP}"

GENERATOR=""
if command -v ninja >/dev/null 2>&1; then
    GENERATOR="-G Ninja"
fi

first=1
for BACKEND in MESHIO NATIVE KRATOS; do
    tree="$SOURCE_DIR/build/bench-$(echo "$BACKEND" | tr '[:upper:]' '[:lower:]')"
    echo "== $BACKEND: configure + build ($tree) =="
    # shellcheck disable=SC2086
    # HDF5/netCDF off: none of the benchmarked formats need them, and it
    # keeps the three throwaway trees small and dependency-free.
    cmake $GENERATOR -S "$SOURCE_DIR" -B "$tree" \
        -DCMAKE_BUILD_TYPE=Release \
        -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
        -DMESHIOPLUSPLUS_BUILD_BENCHMARKS=ON \
        -DMESHIOPLUSPLUS_MESH_BACKEND="$BACKEND" \
        -DMESHIOPLUSPLUS_PARALLEL_BACKEND="$PARALLEL_BACKEND" \
        -DMESHIOPLUSPLUS_WITH_HDF5=OFF \
        -DMESHIOPLUSPLUS_WITH_NETCDF=OFF \
        >/dev/null
    cmake --build "$tree" --target meshioplusplus_bench -j >/dev/null
    echo "== $BACKEND: run (n=$N) =="
    if [ "$first" = 1 ]; then
        "$tree/meshioplusplus_bench" "$N" > "$OUT"
        first=0
    else
        "$tree/meshioplusplus_bench" "$N" | tail -n +2 >> "$OUT"
    fi
done

echo
echo "Results written to $OUT:"
column -s, -t < "$OUT"
