#!/bin/sh
# Run meshioplusplus_bench_containers (roadmap 5.22.1) with one
# (workload, container) per process, so the peak resident set of each row
# belongs to it, and collate one CSV.
#
#   tools/bench_containers.sh BINARY [OUT.csv] [-- bench_containers args]
#
# BINARY is a build with -DMESHIOPLUSPLUS_BUILD_BENCHMARKS=ON. The CSV is
# `workload,container,n,median_s,checksum,maxrss_kb`; the script fails when the
# checksums of one (workload, n) differ between containers.
set -eu
BIN=${1:?usage: bench_containers.sh BINARY [OUT.csv] [-- args]}
OUT=${2:-benchmark/results_containers.csv}
shift $(( $# > 2 ? 2 : $# ))
[ "${1:-}" = "--" ] && shift
"$BIN" --n 1000 --runs 1 "$@" | tail -n +2 | cut -d, -f1,2 | sort -u > "$OUT.pairs"
echo "workload,container,n,median_s,checksum,maxrss_kb" > "$OUT"
while IFS=, read -r workload container; do
    "$BIN" --workload "$workload" --container "$container" "$@" | tail -n +2 >> "$OUT"
done < "$OUT.pairs"
rm -f "$OUT.pairs"
awk -F, 'NR > 1 { k = $1 "," $3; if (k in c && c[k] != $5) { print "checksum mismatch: " $0; bad = 1 } c[k] = $5 } END { exit bad }' "$OUT"
