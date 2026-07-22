#!/bin/sh
# amalgamate.sh -- regenerate the single-header amalgamation of the meshio++ C++
# core at src/single_include/meshioplusplus/meshioplusplus.hpp.
#
#   ./tools/amalgamate.sh              # regenerate the committed single header
#   ./tools/amalgamate.sh --check      # fail if the committed header is stale
#   ./tools/amalgamate.sh --smoke      # regenerate + smoke-compile the header
#
# The single header is committed to the repo (a first-class drop-in deliverable,
# like nlohmann/json). CI runs `--check` (and `--smoke`) so a PR that changes
# src/cpp/ without regenerating the header fails.

set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE_DIR=$(dirname -- "$SCRIPT_DIR")

OUTPUT="$SOURCE_DIR/src/single_include/meshioplusplus/meshioplusplus.hpp"
PYTHON="${PYTHON:-python3}"
DO_CHECK="no"
DO_SMOKE="no"

usage() {
    cat <<EOF
Usage: $0 [options]
  --check     regenerate to a temp file and fail if it differs from the
              committed single header (does not modify the tree)
  --smoke     after regenerating, smoke-compile the header in decls-only and
              MESHIOPLUSPLUS_IMPLEMENTATION modes and link two TUs together
  --output <path>   output header path (default: src/single_include/.../meshioplusplus.hpp)
  -h, --help  this help
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --check) DO_CHECK="yes"; shift ;;
        --smoke) DO_SMOKE="yes"; shift ;;
        --output) OUTPUT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

GEN="$SCRIPT_DIR/amalgamate/amalgamate.py"

if [ "$DO_CHECK" = "yes" ]; then
    TMP=$(mktemp)
    trap 'rm -f "$TMP"' EXIT
    "$PYTHON" "$GEN" --repo-root "$SOURCE_DIR" --output "$TMP"
    if ! diff -u "$OUTPUT" "$TMP" >/dev/null 2>&1; then
        echo "error: $OUTPUT is stale." >&2
        echo "       Run ./tools/amalgamate.sh and commit the result." >&2
        diff -u "$OUTPUT" "$TMP" || true
        exit 1
    fi
    echo "amalgamate: single header is up to date."
    exit 0
fi

"$PYTHON" "$GEN" --repo-root "$SOURCE_DIR" --output "$OUTPUT"
echo "amalgamate: wrote $OUTPUT"

if [ "$DO_SMOKE" = "yes" ]; then
    CXX="${CXX:-g++}"
    INCDIR=$(dirname -- "$(dirname -- "$OUTPUT")")   # the src/single_include/ dir
    PUGI="$SOURCE_DIR/src/cpp/third_party/pugixml"
    WORK=$(mktemp -d)
    trap 'rm -rf "$WORK"' EXIT

    # 1) declarations-only TU (safe to include in many TUs)
    cat > "$WORK/decl_a.cpp" <<EOF
#include "meshioplusplus/meshioplusplus.hpp"
EOF
    # 2) the single implementation TU
    cat > "$WORK/impl.cpp" <<EOF
#define MESHIOPLUSPLUS_IMPLEMENTATION
#include "meshioplusplus/meshioplusplus.hpp"
EOF
    # 3) a second decls-only TU exercising the public API, with a main()
    cat > "$WORK/main.cpp" <<EOF
#include <cstdio>
#include "meshioplusplus/meshioplusplus.hpp"
int main() {
    // Touch the public API surface so the declarations are actually used and
    // linked against the implementation TU.
    const auto& readers = meshioplusplus::registry_readers();
    const std::string fmt = meshioplusplus::resolve_format("mesh.vtk", "");
    std::printf("readers=%zu fmt=%s\n", readers.size(), fmt.c_str());
    return 0;
}
EOF
    echo "amalgamate: smoke-compiling (decls-only, implementation, two-TU link)..."
    STD="-std=c++20"
    INCS="-I$INCDIR -I$PUGI"
    "$CXX" $STD $INCS -fsyntax-only "$WORK/decl_a.cpp"
    "$CXX" $STD $INCS -c "$WORK/impl.cpp" -o "$WORK/impl.o"
    "$CXX" $STD $INCS -c "$WORK/main.cpp" -o "$WORK/main.o"
    "$CXX" $STD "$WORK/impl.o" "$WORK/main.o" -o "$WORK/smoke"
    echo "amalgamate: smoke build OK ($WORK/smoke)"
fi
