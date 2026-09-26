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
#
# With `--hash` among the bench_ops args, every row carries a digest of its
# result, and the script fails unless each (op, cells) row has the same digest
# on every backend and thread count (roadmap §3: the determinism check). With
# BASELINE=<earlier.csv> set, the digests must also equal that file's, which is
# how a change that must not alter output proves it.
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

# The determinism check: one digest per (op, cells) across the whole sweep,
# and, with BASELINE, the same digest as the earlier run. Rows without a
# digest (no --hash) are not checked.
status=0
awk -F, -v base="${BASELINE:-}" '
    FNR == 1 { next }
    FILENAME == base { if ($7 != "") ref[$3 "," $4] = $7; next }
    $7 == "" { next }
    {
        key = $3 "," $4
        if (!(key in seen)) { seen[key] = $7; who[key] = $1 "/" $2 }
        else if (seen[key] != $7) {
            printf "DIGEST MISMATCH %s: %s at %s, %s at %s/%s\n", key, seen[key], who[key], $7, $1, $2
            bad = 1
        }
        if (base != "" && (key in ref) && ref[key] != $7) {
            printf "BASELINE MISMATCH %s: %s now at %s/%s, %s in %s\n", key, $7, $1, $2, ref[key], base
            bad = 1
        }
    }
    END { exit bad }
' ${BASELINE:+"$BASELINE"} "$OUT" >&2 || status=1
[ $status -eq 0 ] && grep -q ',[0-9a-f]\{16\}$' "$OUT" && echo "digests agree" >&2
exit $status
