#!/usr/bin/env bash
# Runs a command with a private TMPDIR, removed afterwards.
#
# ctest starts one process per gtest case, and several helpers name their scratch directory
# from a per-process counter, an unseeded rand() or a time-based seed, which repeat across
# processes. `ctest -j` (used by the instrumented coverage run) then lets two tests share
# and delete one directory. Registering this script as CMAKE_CROSSCOMPILING_EMULATOR gives
# every test process its own temp directory, so no helper can collide.
set -u
d=$(mktemp -d)
TMPDIR="$d" "$@"
rc=$?
rm -rf "$d"
exit "$rc"
