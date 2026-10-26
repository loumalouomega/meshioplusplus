#!/bin/sh
# Syntax-checks the sources that include OpenCASCADE headers against one
# OpenCASCADE include directory, so an API break in a new release is caught
# before the release packages build it.
#
# Usage: sh tools/occt_syntax.sh <opencascade include dir>
#
# Deprecation warnings are printed but do not fail the check: they are
# migrated separately. Errors fail it.
set -eu

occt="$1"
root=$(cd "$(dirname "$0")/.." && pwd)
status=0

for src in \
    src/cpp/src/detail/occt_brep.cpp \
    src/cpp/src/formats/cad_read.cpp \
    src/cpp/src/operations/fit_cad.cpp
do
    echo "== $src"
    g++ -std=c++20 -fsyntax-only -DMESHIOPLUSPLUS_HAS_OCCT \
        -isystem "$occt" \
        -I "$root/src/cpp/include" \
        -I "$root/src/cpp/third_party/json/single_include" \
        -I "$root/src/cpp/third_party/pugixml" \
        "$root/$src" || status=1
done

exit $status
