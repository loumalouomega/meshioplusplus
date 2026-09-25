#!/bin/sh
# Build meshioplusplus_bench_ops once per parallel backend and sweep thread
# counts (OMP_NUM_THREADS, which bench_ops also applies to TBB), collating one CSV (doc/benchmarks.md, roadmap §3).
#
#   tools/bench_ops.sh [OUT.csv] [BACKENDS] [THREADS] [-- bench_ops args]
#
# BACKENDS defaults to "SEQ OPENMP TBB" (a backend whose configure fails, e.g.
# TBB not installed, is skipped with a note); THREADS to "1 2 4 8". The SEQ
# backend runs once. Trees go to build/bench-ops-<backend>; set
# CMAKE_BUILD_PARALLEL_LEVEL to cap the build.
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(dirname -- "$HERE")
OUT=${1:-benchmark/results_ops.csv}
BACKENDS=${2:-"SEQ OPENMP TBB"}
THREADS=${3:-"1 2 4 8"}
shift $(( $# > 3 ? 3 : $# ))
[ "${1:-}" = "--" ] && shift
: > "$OUT.tmp"
header=""
for b in $BACKENDS; do
    tree="$ROOT/build/bench-ops-$(echo "$b" | tr '[:upper:]' '[:lower:]')"
    if ! cmake -S "$ROOT" -B "$tree" -G Ninja -DCMAKE_BUILD_TYPE=Release \
        -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF -DMESHIOPLUSPLUS_BUILD_BENCHMARKS=ON \
        -DMESHIOPLUSPLUS_PARALLEL_BACKEND="$b" -DMESHIOPLUSPLUS_WITH_HDF5=OFF \
        -DMESHIOPLUSPLUS_WITH_NETCDF=OFF > "$tree.configure.log" 2>&1; then
        echo "skip $b: configure failed (see $tree.configure.log)" >&2
        continue
    fi
    cmake --build "$tree" --target meshioplusplus_bench_ops
    sweep=$THREADS
    [ "$b" = SEQ ] && sweep=1
    for t in $sweep; do
        echo "== $b, $t thread(s)" >&2
        OMP_NUM_THREADS=$t "$tree/meshioplusplus_bench_ops" "$@" > "$OUT.run"
        if [ -z "$header" ]; then
            header=$(head -1 "$OUT.run")
            echo "$header" >> "$OUT.tmp"
        fi
        tail -n +2 "$OUT.run" >> "$OUT.tmp"
    done
done
rm -f "$OUT.run"
mv "$OUT.tmp" "$OUT"
echo "wrote $OUT" >&2
