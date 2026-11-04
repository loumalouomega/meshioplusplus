#!/bin/sh
# Development build of the MEX gateway into a private directory.
#
#   tools/dev_build.sh OUTDIR [FILE.c ...]
#
# Compiles the always-present core (util.c, mesh.c, mio_mex.c) plus the named
# src/ files into OUTDIR/meshioplusplus_mex.mex. Every command table the build
# does NOT include gets an empty generated stub, so a half-finished port in
# another file cannot break this build. With no FILE arguments, every src/*.c
# is compiled (the full build). The C library is found through
# MESHIOPLUSPLUS_HOME or pkg-config, as in mio_build.m.
set -e
here=$(cd "$(dirname "$0")/.." && pwd)
out=$1
shift || true
[ -n "$out" ] || { echo "usage: $0 OUTDIR [FILE.c ...]" >&2; exit 2; }
mkdir -p "$out"
if [ -n "$MESHIOPLUSPLUS_HOME" ]; then
    cflags="-I$MESHIOPLUSPLUS_HOME/include"
    libs="-L$MESHIOPLUSPLUS_HOME/lib -lmeshioplusplus"
else
    cflags=$(pkg-config --cflags meshioplusplus)
    libs=$(pkg-config --libs meshioplusplus)
fi
tables="mesh ops_core ops_combine ops_data format_info sequence xdmf_series exodus_series render"
if [ $# -eq 0 ]; then
    srcs="$here/src/*.c"
else
    srcs="$here/src/util.c $here/src/mesh.c $here/src/mio_mex.c"
    for f in "$@"; do srcs="$srcs $here/src/$(basename "$f")"; done
    for t in $tables; do
        case " $srcs " in *"/src/$t.c "*) continue ;; esac
        printf '#include "mio_oct.h"\nconst mio_o_entry *mio_o_table_%s(size_t *n) { *n = 0; return NULL; }\n' "$t" \
            > "$out/stub_$t.c"
        srcs="$srcs $out/stub_$t.c"
    done
fi
# shellcheck disable=SC2086
mkoctfile --mex -Wall -Wextra -I"$here/src" $cflags $srcs $libs -o "$out/meshioplusplus_mex.mex"
echo "built $out/meshioplusplus_mex.mex"
