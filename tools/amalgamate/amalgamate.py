#!/usr/bin/env python3
# amalgamate.py -- generate a single-header, header-only amalgamation of the
# meshio++ C++ core.
#
# This implements the classic "amalgamate" approach (in the spirit of
# edlund/amalgamate and the STB single-file-library convention), producing an
# STB-style single header:
#
#   * declarations are always visible on a plain ``#include``;
#   * the out-of-line implementations of every ``src/cpp/src/*.cpp`` (and the
#     bundled pugixml) are compiled only in the one translation unit that first
#     ``#define``s ``MESHIOPLUSPLUS_IMPLEMENTATION`` before including the header.
#
# Rather than naively inlining at each ``#include`` site, headers are emitted
# ONCE, at top level, in dependency (topological) order, with project
# ``#include`` directives stripped. This is what keeps the amalgamation correct:
# a header that includes a common header from *inside* a preprocessor guard
# (e.g. detail/hdf5_util.hpp includes exceptions.hpp inside
# ``#ifdef MESHIOPLUSPLUS_HAS_HDF5``) can never trap that common header inside
# the guard -- every header stands on its own at file scope.
#
# The mesh backend is the one exception: the three ``backends/*`` headers are
# mutually exclusive (selected by mesh.hpp's ``#if/#elif/#else``), so they are
# inlined at their conditional sites inside mesh.hpp -- only the selected
# backend compiles, exactly as in the normal build.
#
# Optional dependencies (HDF5/netCDF/zlib/Eigen) stay behind their existing
# ``MESHIOPLUSPLUS_HAS_*`` guards; the header compiles with zero external
# dependencies unless the consumer opts in. Output is deterministic (sorted
# traversal) so it is byte-stable, committable, and verifiable in CI.

import argparse
import re
import sys
from pathlib import Path

_QUOTED_INCLUDE = re.compile(r'^\s*#\s*include\s*"([^"]+)"')
_PRAGMA_ONCE = re.compile(r"^\s*#\s*pragma\s+once\b")

_REPO_ROOT = None


def _display(path):
    try:
        return str(path.relative_to(_REPO_ROOT))
    except ValueError:
        return path.name


def strip_leading_banner(text):
    """Drop the leading license banner / ``#pragma once`` from one file.

    Every file opens with an ASCII-art license banner (a run of ``//`` comment
    lines) followed, in headers, by ``#pragma once``. Keeping 60+ copies would
    bloat the header, so drop the leading run of blank/``//``-comment lines plus
    a trailing ``#pragma once``. Block comments and code are preserved.
    """
    lines = text.splitlines(keepends=True)
    i, n = 0, len(lines)
    while i < n:
        s = lines[i].strip()
        if s == "" or s.startswith("//") or _PRAGMA_ONCE.match(lines[i]):
            i += 1
            continue
        break
    return lines[i:]


class Amalgamator:
    def __init__(self, include_paths, inline_files):
        # Directories searched (after the including file's own dir) to resolve a
        # quoted include to a file on disk.
        self.include_paths = [Path(p).resolve() for p in include_paths]
        # Project headers that must be inlined at their include site rather than
        # emitted at top level (the mutually-exclusive mesh backends).
        self.inline_files = {Path(p).resolve() for p in inline_files}
        self.emitted = set()  # absolute paths already written out
        self.out = []

    def resolve(self, quoted, from_dir):
        cand = (from_dir / quoted).resolve()
        if cand.is_file():
            return cand
        for base in self.include_paths:
            cand = (base / quoted).resolve()
            if cand.is_file():
                return cand
        return None

    def project_deps(self, path):
        """Resolved project/pugixml headers this file #includes (any nesting)."""
        deps = []
        from_dir = path.parent
        for line in path.read_text().splitlines():
            m = _QUOTED_INCLUDE.match(line)
            if m:
                t = self.resolve(m.group(1), from_dir)
                if t is not None:
                    deps.append(t)
        return deps

    def emit_body(self, path):
        """Write ``path``'s body: banner + #pragma once + project includes
        stripped; a project include to an *inline* file is expanded in place;
        everything else (system includes, code, this file's own #ifdefs) is kept
        verbatim so each file's preprocessor/namespace structure stays balanced."""
        from_dir = path.parent
        for line in strip_leading_banner(path.read_text()):
            m = _QUOTED_INCLUDE.match(line)
            if m:
                t = self.resolve(m.group(1), from_dir)
                if t is not None:
                    if t in self.inline_files:
                        self.inline_at_site(t)
                    # else: emitted at top level -> drop the include line
                    continue
            if _PRAGMA_ONCE.match(line):
                continue
            self.out.append(line)

    def inline_at_site(self, path):
        path = path.resolve()
        if path in self.emitted:
            return
        self.emitted.add(path)
        self.out.append(f"// ===== begin {_display(path)} =====\n")
        self.emit_body(path)
        self.out.append(f"// ===== end {_display(path)} =====\n")

    def emit_toplevel(self, path):
        """Emit a header once at file scope (used for the topological order)."""
        self.inline_at_site(path)


def topo_order(files, deps_fn):
    """Deterministic post-order (dependencies first) over ``files``."""
    files = list(files)
    fileset = set(files)
    order, visiting, done = [], set(), set()

    def visit(f):
        if f in done:
            return
        if f in visiting:
            return  # cycle: break it (headers are #pragma once; treat as DAG)
        visiting.add(f)
        for d in deps_fn(f):
            if d in fileset:
                visit(d)
        visiting.discard(f)
        done.add(f)
        order.append(f)

    for f in sorted(files):
        visit(f)
    return order


def main(argv=None):
    global _REPO_ROOT
    ap = argparse.ArgumentParser(
        description="Amalgamate the meshio++ C++ core into one header."
    )
    ap.add_argument("--repo-root", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args(argv)

    repo = Path(args.repo_root).resolve()
    _REPO_ROOT = repo
    include_dir = repo / "src" / "cpp" / "include"
    pugixml_dir = repo / "src" / "cpp" / "third_party" / "pugixml"

    all_headers = sorted((include_dir / "meshioplusplus").rglob("*.hpp"))
    # Every header is emitted once, at top level, in dependency order. The three
    # backend structs have distinct names (Mesh/CellBlock, NativeMesh/..., and
    # KratosMesh + the ModelPart classes), so they coexist; mesh.hpp's
    # #if/#elif/#else then only selects which one `Mesh` aliases. Emitting them
    # unconditionally (rather than trapping them inside mesh.hpp's branches) is
    # also what lets the always-emitted kratos_bridge.hpp see ModelPart.
    amalg = Amalgamator(include_paths=[include_dir, pugixml_dir], inline_files=set())

    order = topo_order(all_headers, amalg.project_deps)

    amalg.out.append(BANNER)
    amalg.out.append("#pragma once\n\n")
    amalg.out.append(DEFAULTS)
    amalg.out.append("\n")

    amalg.out.append(
        "// ================= DECLARATIONS (always compiled) =================\n"
    )
    for h in order:
        amalg.emit_toplevel(h)

    # Implementation section: bundled pugixml + every core .cpp. Any third-party
    # header the sources need but the decl section didn't emit (pugixml.hpp,
    # pugiconfig.hpp) is emitted here first, at top level.
    sources = [pugixml_dir / "pugixml.cpp"]
    sources += sorted((repo / "src" / "cpp" / "src").rglob("*.cpp"))

    prelude = []
    for s in sources:
        for d in amalg.project_deps(s):
            if (
                d not in amalg.emitted
                and d not in prelude
                and d.suffix in (".hpp", ".h")
            ):
                prelude.append(d)
    prelude = topo_order(prelude, amalg.project_deps)

    amalg.out.append("\n#ifdef MESHIOPLUSPLUS_IMPLEMENTATION\n")
    amalg.out.append("// ================= IMPLEMENTATION =================\n")
    for h in prelude:
        amalg.emit_toplevel(h)
    for s in sources:
        amalg.emit_toplevel(s)
    amalg.out.append("#endif  // MESHIOPLUSPLUS_IMPLEMENTATION\n")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(amalg.out))
    print(
        f"amalgamate: {len(all_headers)} headers + {len(sources)} sources "
        f"-> {out_path} ({out_path.stat().st_size} bytes)",
        file=sys.stderr,
    )


BANNER = """\
// meshio++ -- single-header, header-only amalgamation of the C++ core.
//
//  License: MIT (meshio++ default license: LICENSE). Bundles pugixml (MIT).
//  Main authors: Vicente Mataix Ferrandiz
//
//  *** GENERATED FILE -- DO NOT EDIT BY HAND. ***
//  Regenerate with:  ./tools/amalgamate.sh
//  (CI verifies this file is up to date; edit the sources under src/cpp/, not here.)
//
//  Usage (STB-style, header-only):
//
//      // in exactly ONE translation unit:
//      #define MESHIOPLUSPLUS_IMPLEMENTATION
//      #include "meshioplusplus/meshioplusplus.hpp"
//
//      // in every other translation unit -- declarations only:
//      #include "meshioplusplus/meshioplusplus.hpp"
//
//  Mesh backend defaults to MESHIO, parallel backend to sequential. Optional
//  formats stay off unless you define the matching macro AND link the library:
//      MESHIOPLUSPLUS_HAS_HDF5  (CGNS/HMF/H5M/MED/XDMF-HDF)  -> link hdf5
//      MESHIOPLUSPLUS_HAS_NETCDF (Exodus)                    -> link netcdf
//      MESHIOPLUSPLUS_HAS_ZLIB  (VTU zlib compression)       -> link z
//      MESHIOPLUSPLUS_HAS_EIGEN (MED transpose fast path)    -> add Eigen to the include path
//      MESHIOPLUSPLUS_HAS_KAHIP (partition kahip backend)    -> link kahip
"""

DEFAULTS = """\
#if !defined(MESHIOPLUSPLUS_PARALLEL_SEQ) && !defined(MESHIOPLUSPLUS_PARALLEL_STL) && \\
    !defined(MESHIOPLUSPLUS_PARALLEL_OPENMP) && !defined(MESHIOPLUSPLUS_PARALLEL_TBB)
#define MESHIOPLUSPLUS_PARALLEL_SEQ
#endif
// Mesh backend: MESHIO is the no-macro default (see mesh.hpp's #else); define
// MESHIOPLUSPLUS_MESH_BACKEND_NATIVE or _KRATOS before including to change it.
"""


if __name__ == "__main__":
    main()
