#!/usr/bin/env python3
"""Convert an lcov tracefile to SonarQube's generic test-coverage XML.

The ``coverage`` job in ``.github/workflows/ci.yml`` captures the C++ core's
coverage once with lcov (``cpp-coverage.info``, already filtered of third-party
and unreachable backend code) and uploads it to Codecov. SonarCloud's C/C++
analyzer does not read lcov, so this rewrites the same tracefile in Sonar's
generic format (``sonar.coverageReportPaths``) instead of re-running gcov: both
services then see the same filtered denominator (doc/static_analysis.md).

    tools/lcov_to_sonar.py cpp-coverage.info cpp-coverage-sonar.xml --root .

Paths are written relative to ``--root``; files outside it (generated sources
in a build tree, system headers) are dropped, since Sonar only maps coverage
onto files it indexes. Records of the same file are merged: line hits are
summed and a branch counts as covered when any record took it. Standard
library only.
"""

import argparse
import os
import sys
from xml.sax.saxutils import quoteattr


def parse_lcov(lines):
    """Return ``{path: (line_hits, branches)}`` from lcov tracefile lines.

    ``line_hits`` maps a line number to its summed execution count;
    ``branches`` maps a line number to ``{(block, branch): taken}``, where
    ``taken`` is the summed count (``-``, a branch whose block never ran,
    counts as 0).
    """
    files = {}
    current = None
    for raw in lines:
        line = raw.strip()
        if line.startswith("SF:"):
            path = line[3:]
            current = files.setdefault(path, ({}, {}))
        elif line == "end_of_record":
            current = None
        elif current is None:
            continue
        elif line.startswith("DA:"):
            fields = line[3:].split(",")
            lineno, hits = int(fields[0]), int(fields[1])
            current[0][lineno] = current[0].get(lineno, 0) + hits
        elif line.startswith("BRDA:"):
            lineno, block, branch, taken = line[5:].split(",")[:4]
            count = 0 if taken == "-" else int(taken)
            per_line = current[1].setdefault(int(lineno), {})
            key = (block, branch)
            per_line[key] = per_line.get(key, 0) + count
    return files


def to_sonar_xml(files, root):
    """Render parsed lcov data as Sonar generic coverage XML (a string)."""
    root = os.path.abspath(root)
    out = ['<?xml version="1.0" encoding="UTF-8"?>', '<coverage version="1">']
    for path in sorted(files):
        absolute = os.path.abspath(os.path.join(root, path))
        relative = os.path.relpath(absolute, root)
        if relative == os.pardir or relative.startswith(os.pardir + os.sep):
            continue
        line_hits, branches = files[path]
        lines = sorted(set(line_hits) | set(branches))
        if not lines:
            continue
        out.append(f"  <file path={quoteattr(relative.replace(os.sep, '/'))}>")
        for lineno in lines:
            if lineno <= 0:
                continue
            hits = line_hits.get(lineno, 0)
            attrs = f'lineNumber="{lineno}" covered="{"true" if hits > 0 else "false"}"'
            if lineno in branches:
                taken = branches[lineno].values()
                covered = sum(1 for t in taken if t > 0)
                attrs += f' branchesToCover="{len(taken)}" coveredBranches="{covered}"'
            out.append(f"    <lineToCover {attrs}/>")
        out.append("  </file>")
    out.append("</coverage>")
    return "\n".join(out) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("tracefile", help="lcov tracefile (.info) to read")
    parser.add_argument("output", help="Sonar generic coverage XML to write")
    parser.add_argument(
        "--root", default=".", help="repository root the written paths are relative to"
    )
    args = parser.parse_args(argv)
    with open(args.tracefile, encoding="utf-8", errors="replace") as fh:
        files = parse_lcov(fh)
    xml = to_sonar_xml(files, args.root)
    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(xml)
    kept = xml.count("<file ")
    print(f"wrote {args.output}: {kept} of {len(files)} files", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
