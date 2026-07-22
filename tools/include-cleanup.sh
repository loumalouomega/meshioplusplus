#!/bin/sh
# include-cleanup.sh -- IWYU-style include hygiene for the meshio++ C++ core,
# driven by clang-tidy's misc-include-cleaner check (see .clang-tidy).
#
#   ./tools/include-cleanup.sh --check   # CI gate: fail on any UNUSED include
#   ./tools/include-cleanup.sh --fix     # remove unused includes from the tree
#
# Policy (deliberately conservative -- "include only what is really necessary"):
#   * UNUSED includes (a header whose symbols this TU never uses) are the gate:
#     --check fails on them, --fix removes them.
#   * MISSING includes (a symbol used but only pulled in transitively) are
#     reported as ADVISORY only -- the project's convention is that a format's
#     .cpp may lean on its own format header to pull in the mesh types, so these
#     are informational, never fatal, and never auto-applied.
#
# It needs a CMake compile database (compile_commands.json). By default it
# configures a throwaway build with the optional deps ON so every #ifdef-guarded
# translation unit has a compile command and is analyzed; pass --build-dir to
# reuse an existing configured tree. NOTE: run the gate where HDF5/netCDF are
# installed (as CI does) -- with those libraries absent, the HDF5/netCDF format
# registrations compile out and their headers look spuriously unused.
#
# Only our own sources are analyzed (src/cpp/src, src/cpp/include, bindings/c) -- never
# src/cpp/third_party. Deliberately-kept includes carry `// IWYU pragma: keep`.

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_DIR=$(dirname -- "$SCRIPT_DIR")

CLANG_TIDY="${CLANG_TIDY:-clang-tidy}"
BUILD_DIR=""
MODE="check"

usage() {
    cat <<EOF
Usage: $0 [--fix | --check] [--build-dir <dir>]
  --check             fail (exit 1) if any UNUSED include is found (default; CI)
  --fix               remove unused includes from the tree
  --build-dir <dir>   reuse an existing CMake build with compile_commands.json
                      (default: configure build/include-cleanup-cc)
  -h, --help          this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --fix) MODE="fix"; shift ;;
        --check) MODE="check"; shift ;;
        --build-dir) BUILD_DIR="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [ -z "$BUILD_DIR" ]; then
    BUILD_DIR="$SOURCE_DIR/build/include-cleanup-cc"
    echo "include-cleanup: configuring compile database in $BUILD_DIR ..."
    # BUILD_PYTHON=OFF: the pybind boundary (bindings/python/np_conversions.hpp)
    # is the sanctioned uniform-API exception and is not analyzed here. C_API=ON
    # pulls bindings/c into the database. SEQ backend avoids needing <omp.h>/TBB
    # on the clang-tidy include path (include analysis is identical either way).
    cmake -S "$SOURCE_DIR" -B "$BUILD_DIR" -G Ninja \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
        -DMESHIOPLUSPLUS_BUILD_C_API=ON \
        -DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ >/dev/null
fi

if [ ! -f "$BUILD_DIR/compile_commands.json" ]; then
    echo "error: no compile_commands.json in $BUILD_DIR" >&2
    exit 1
fi

LOG=$(mktemp)
trap 'rm -f "$LOG"' EXIT

# Analyze every core translation unit one at a time (a single clang-tidy call
# over many files trips a "no input files" quirk); include-cleaner also visits
# the headers each TU pulls in (filtered by HeaderFilterRegex in .clang-tidy).
n=0
for f in $(find "$SOURCE_DIR/src/cpp/src" "$SOURCE_DIR/bindings/c" -name '*.cpp' | sort); do
    "$CLANG_TIDY" -p "$BUILD_DIR" --quiet "$f" 2>&1 >/dev/null | tee -a "$LOG" >/dev/null || true
    n=$((n + 1))
done
echo "include-cleanup: analyzed $n translation units."

python3 - "$LOG" "$MODE" <<'PY'
import re, sys
log_path, mode = sys.argv[1], sys.argv[2]
unused, adds = [], 0
for line in open(log_path):
    m = re.search(r'(\S+):(\d+):\d+: warning: included header (\S+) is not used', line)
    if m:
        unused.append((m.group(1), int(m.group(2)), m.group(3)))
        continue
    if 'no header providing' in line:
        adds += 1

if mode == "fix":
    by_file = {}
    for path, ln, hdr in unused:
        by_file.setdefault(path, []).append((ln, hdr))
    total = 0
    for path, items in by_file.items():
        lines = open(path).read().splitlines(keepends=True)
        for ln, hdr in sorted(items, reverse=True):
            if hdr.split('/')[-1] in lines[ln - 1]:
                del lines[ln - 1]
                total += 1
        open(path, "w").writelines(lines)
    print(f"include-cleanup: removed {total} unused include(s).")
    if adds:
        print(f"include-cleanup: {adds} missing-include suggestion(s) NOT applied "
              f"(advisory; add '// IWYU pragma: keep' or the include by hand if wanted).")
    sys.exit(0)

# check mode
if unused:
    print(f"include-cleanup: {len(unused)} UNUSED include(s) found:")
    for path, ln, hdr in unused:
        print(f"  {path}:{ln}: unused <{hdr}>")
    print("Run ./tools/include-cleanup.sh --fix to remove them.")
    sys.exit(1)
print("include-cleanup: no unused includes.")
if adds:
    print(f"include-cleanup: {adds} missing-include suggestion(s) (advisory only).")
sys.exit(0)
PY
