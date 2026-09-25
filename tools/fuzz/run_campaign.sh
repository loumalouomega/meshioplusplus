#!/bin/bash
# Run the reader fuzz target over every registry format (doc/fuzzing.md).
#
#   tools/fuzz/run_campaign.sh BUILD_DIR SEEDS_DIR OUT_DIR [SECONDS_PER_FORMAT]
#
# BUILD_DIR holds meshioplusplus_fuzz_read and meshioplusplus_fuzz_replay (a
# --fuzzers tree); SEEDS_DIR is tools/fuzz/seed_corpus.py's output. Each format
# runs for SECONDS_PER_FORMAT (default 60) with its dictionary when one exists;
# findings land in OUT_DIR/<format>/ (crash-*, leak-*, timeout-*, oom-*) and a
# one-line verdict per format in OUT_DIR/summary.txt. FORMATS (space-separated)
# narrows the list; SHARD/NSHARDS (0-based) split it for a CI matrix.
# Exit status: 1 when any format produced a finding.
set -u
BUILD=$1
SEEDS=$2
OUT=$3
SECS=${4:-60}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
SKIP=$(grep -v '^#' "$HERE/not_fuzzed.txt" | grep -v '^$')
ALL=${FORMATS:-$("$BUILD/meshioplusplus_fuzz_replay" -list-formats)}
SHARD=${SHARD:-0}
NSHARDS=${NSHARDS:-1}
: > "$OUT/summary.txt"
status=0
i=0
for fmt in $ALL; do
    idx=$i
    i=$((i + 1))
    [ $((idx % NSHARDS)) -eq "$SHARD" ] || continue
    if echo "$SKIP" | grep -qx "$fmt"; then
        continue
    fi
    dir="$OUT/$fmt"
    mkdir -p "$dir/corpus"
    [ -d "$SEEDS/$fmt" ] && cp -n "$SEEDS/$fmt"/* "$dir/corpus/" 2>/dev/null
    dict=()
    [ -f "$HERE/dicts/$fmt.dict" ] && dict=(-dict="$HERE/dicts/$fmt.dict")
    MIO_FUZZ_FORMAT=$fmt "$BUILD/meshioplusplus_fuzz_read" "${dict[@]}" \
        -max_total_time="$SECS" -rss_limit_mb=2048 -malloc_limit_mb=2048 -timeout=10 \
        -max_len=65536 -print_final_stats=1 -artifact_prefix="$dir/" \
        "$dir/corpus" > "$dir/fuzz.log" 2>&1
    rc=$?
    found=$(ls "$dir" | grep -E '^(crash|leak|timeout|oom|slow-unit)-' | head -1)
    if [ -n "$found" ] || [ $rc -ne 0 ]; then
        echo "$fmt: FINDING ${found:-exit $rc}" | tee -a "$OUT/summary.txt"
        status=1
    else
        execs=$(grep -ao 'stat::number_of_executed_units: *[0-9]*' "$dir/fuzz.log" | grep -o '[0-9]*$')
        echo "$fmt: ok (${execs:-?} execs)" | tee -a "$OUT/summary.txt"
    fi
done
exit $status
