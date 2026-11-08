#!/usr/bin/env python3
"""Convert an lcov tracefile to SonarQube's generic coverage XML (doc/static_analysis.md).

The ``coverage`` job of ``ci.yml`` captures the C++ coverage with lcov and
uploads it to Codecov. SonarCloud's C/C++ analyzer does not read lcov, but every
SonarQube analyzer reads the generic format, so this turns the same, already
filtered tracefile into it. Both reports then share one denominator (the
``lcov --remove`` list of that job).

    tools/lcov_to_sonar.py cpp-coverage.info cpp-coverage-sonar.xml --root "$PWD"

A source path under ``--root`` is written relative to it, which is how the
scanner names files; a source outside it (a system header lcov kept) is left
out. ``DA`` records give the lines and their hits, ``BRDA`` records the branches
of a line. Standard library only.
"""

import argparse
import os
import sys
from xml.sax.saxutils import quoteattr


def parse_lcov(lines):
    """``{path: {line: [hits, branches, covered_branches]}}`` of a tracefile."""
    files = {}
    current = None
    for raw in lines:
        line = raw.strip()
        if line.startswith("SF:"):
            current = files.setdefault(line[3:], {})
        elif line == "end_of_record":
            current = None
        elif current is None:
            continue
        elif line.startswith("DA:"):
            fields = line[3:].split(",")
            number, hits = int(fields[0]), int(float(fields[1]))
            entry = current.setdefault(number, [0, 0, 0])
            entry[0] += max(hits, 0)
        elif line.startswith("BRDA:"):
            number, _block, _branch, taken = line[5:].split(",", 3)
            entry = current.setdefault(int(number), [0, 0, 0])
            entry[1] += 1
            if taken not in ("-", "0"):
                entry[2] += 1
    return files


def relative_to(path, root):
    """``path`` relative to ``root`` with ``/`` separators, or ``None`` outside it."""
    full = os.path.normpath(path if os.path.isabs(path) else os.path.join(root, path))
    rel = os.path.relpath(full, root)
    if rel == os.pardir or rel.startswith(os.pardir + os.sep):
        return None
    return rel.replace(os.sep, "/")


def to_generic_xml(files, root):
    """The generic coverage document of parsed lcov ``files``."""
    out = ['<coverage version="1">']
    for path in sorted(files):
        rel = relative_to(path, root)
        if rel is None or not files[path]:
            continue
        out.append(f"  <file path={quoteattr(rel)}>")
        for number in sorted(files[path]):
            hits, branches, covered = files[path][number]
            attrs = f'lineNumber="{number}" covered="{"true" if hits > 0 else "false"}"'
            if branches:
                attrs += f' branchesToCover="{branches}" coveredBranches="{covered}"'
            out.append(f"    <lineToCover {attrs}/>")
        out.append("  </file>")
    out.append("</coverage>")
    return "\n".join(out) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("tracefile", help="the lcov .info file")
    parser.add_argument("output", help="the generic coverage XML to write")
    parser.add_argument(
        "--root", default=os.getcwd(), help="the repository root (default: cwd)"
    )
    args = parser.parse_args(argv)
    with open(args.tracefile, encoding="utf-8", errors="replace") as fh:
        files = parse_lcov(fh)
    root = os.path.abspath(args.root)
    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(to_generic_xml(files, root))
    kept = sum(1 for p in files if relative_to(p, root) is not None and files[p])
    print(
        f"lcov_to_sonar: {kept} of {len(files)} source files written to {args.output}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
