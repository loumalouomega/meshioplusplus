#!/usr/bin/env python3
"""A local code-duplication report, CPD-style (doc/static_analysis.md).

SonarCloud measures duplication on ``master`` only. This gives the same kind of
number on any checkout, so a cleaning change can show its effect before it
lands: every source is tokenised (C++: comments and ``#include`` lines dropped;
Python: through ``tokenize``), identifiers, numbers and strings are normalised
so a renamed copy still matches, and every window of ``--tokens`` consecutive
tokens is hashed. A window seen in two places marks its lines as duplicated.

    tools/dup_report.py                    # C++ and Python, top 25 of each list
    tools/dup_report.py --lang cpp --top 60
    tools/dup_report.py --json report.json

The report gives the duplicated-line density per language, the files with the
most duplicated lines, and the largest regions shared by two files with their
line ranges. Excluded: vendored code, the amalgamated single header and the
generated tables (the list ``sonar-project.properties`` excludes), and the
files that are tables alone (its ``sonar.cpd.exclusions``). Windows repeated
more than ``--max-occurrences`` times (long runs of identical table rows) are
ignored, and so is a window made only of literals and table punctuation: once
strings and numbers are normalised, two literal tables of the same shape match
although no entry agrees (roadmap §8). The R binding's C sources are scanned
with the C++ ones, as SonarCloud scans them. Standard library only; it reads
files and writes nothing unless ``--json`` is given.
"""

import argparse
import collections
import io
import json
import keyword
import os
import re
import sys
import tokenize

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CPP_ROOTS = (
    "src/cpp/src",
    "src/cpp/include/meshioplusplus",
    "src/cpp/cli",
    "src/cpp/providers",
    "bindings/c",
    "bindings/python",
    "bindings/r/meshioplusplus/src",
)
PY_ROOTS = ("src/python/meshioplusplus",)
EXCLUDED_PARTS = ("third_party", "single_include", "node_modules")
EXCLUDED_FILES = {
    "hex_dual_scheme_tables.cpp",
    "colormap.cpp",
    "colormap.hpp",
    "_colormap.py",
    "_format_table.py",
    "_schema.py",
    # Tables alone (sonar.cpd.exclusions): per-format node permutations.
    "node_order.cpp",
}

_CPP_TOKEN = re.compile(
    r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[A-Za-z_]\w*|\d[\w.]*|::|->|<<|>>'
    r"|[-+*/%=!<>&|^]=?|&&|\|\||[{}()\[\];,.?:~#]"
)
# Kept verbatim so that structure, not just shape, has to match.
_CPP_KEYWORDS = frozenset(
    "if else for while do return switch case default break continue const auto "
    "static inline constexpr struct class namespace using template typename "
    "public private protected void double int float bool char true false "
    "nullptr new delete std size_t int64_t uint8_t int32_t vector array".split()
)


# What a literal table is made of once strings and numbers are normalised: a
# window holding almost nothing else is a run of table rows, not copied code
# (roadmap §8). One token in ten may be something else, so the windows that
# straddle a table's declaration (``const T kName[] = {``) go with its rows.
# The boolean and null literals of both languages count as literals: an option
# table's ``{"name", {}, true},`` rows are a table, not code.
_TABLE_TOKENS = frozenset(
    ("N", "S", "{", "}", "[", "]", "(", ")", ",", ":", ";", "-", "+", ".")
    + ("true", "false", "nullptr", "True", "False", "None")
)


def is_table_window(tokens):
    """True when ``tokens`` is literals and table punctuation, bar one in ten."""
    other = sum(1 for t in tokens if t not in _TABLE_TOKENS)
    return other <= len(tokens) // 10


def cpp_tokens(text):
    """``[(normalised token, line)]`` of a C/C++ source."""
    text = re.sub(
        r"/\*.*?\*/", lambda m: "\n" * m.group(0).count("\n"), text, flags=re.S
    )
    out = []
    for lineno, line in enumerate(text.split("\n"), 1):
        line = re.sub(r"//.*", "", line)
        if line.lstrip().startswith("#include"):
            continue
        for tok in _CPP_TOKEN.findall(line):
            if tok[0].isalpha() or tok[0] == "_":
                tok = tok if tok in _CPP_KEYWORDS else "I"
            elif tok[0].isdigit():
                tok = "N"
            elif tok[0] in "\"'":
                tok = "S"
            out.append((tok, lineno))
    return out


_PY_SKIP = {
    tokenize.COMMENT,
    tokenize.NL,
    tokenize.NEWLINE,
    tokenize.INDENT,
    tokenize.DEDENT,
    tokenize.ENCODING,
    tokenize.ENDMARKER,
}


def py_tokens(text):
    """``[(normalised token, line)]`` of a Python source (empty if it does not tokenise)."""
    out = []
    try:
        for tok in tokenize.generate_tokens(io.StringIO(text).readline):
            if tok.type in _PY_SKIP:
                continue
            if tok.type == tokenize.NAME:
                norm = tok.string if keyword.iskeyword(tok.string) else "I"
            elif tok.type == tokenize.NUMBER:
                norm = "N"
            elif tok.type == tokenize.STRING:
                norm = "S"
            else:
                norm = tok.string
            out.append((norm, tok.start[0]))
    except (tokenize.TokenError, SyntaxError):
        return []
    return out


def collect(root, roots, suffixes):
    files = []
    for rel in roots:
        for dirpath, dirnames, filenames in os.walk(os.path.join(root, rel)):
            dirnames[:] = sorted(d for d in dirnames if d not in EXCLUDED_PARTS)
            for name in sorted(filenames):
                if name.endswith(suffixes) and name not in EXCLUDED_FILES:
                    files.append(os.path.join(dirpath, name))
    return files


def analyse(paths, tokenizer, window, max_occurrences, root):
    """Duplicated lines per file and the largest shared regions between files."""
    toks = {}
    for path in paths:
        with open(path, encoding="utf-8", errors="replace") as fh:
            toks[os.path.relpath(path, root)] = tokenizer(fh.read())
    index = collections.defaultdict(list)
    for name, seq in toks.items():
        norm = [t for t, _ in seq]
        for i in range(len(norm) - window + 1):
            win = tuple(norm[i : i + window])
            if is_table_window(win):
                continue
            index[hash(win)].append((name, i))
    dup = collections.defaultdict(set)
    pairs = collections.defaultdict(list)
    for occ in index.values():
        if len(occ) < 2 or len(occ) > max_occurrences:
            continue
        for name, i in occ:
            seq = toks[name]
            dup[name].update(seq[k][1] for k in range(i, i + window))
        for a in range(len(occ)):
            for b in range(a + 1, len(occ)):
                (fa, ia), (fb, ib) = occ[a], occ[b]
                if fa == fb:
                    continue
                if fa > fb:
                    fa, ia, fb, ib = fb, ib, fa, ia
                pairs[(fa, fb)].append((ia, ib))
    # Merge consecutive windows on one diagonal into a region.
    regions = []
    for (fa, fb), hits in pairs.items():
        hits.sort(key=lambda h: (h[0] - h[1], h[0]))
        run = None
        for ia, ib in hits:
            if run and ia - ib == run[0] - run[1] and ia <= run[2] + 1:
                run[2] = ia
                continue
            if run:
                regions.append((fa, fb, *run))
            run = [ia, ib, ia]
        if run:
            regions.append((fa, fb, *run))
    out_regions = []
    for fa, fb, ia, ib, ie in regions:
        span = ie - ia
        out_regions.append(
            {
                "tokens": span + window,
                "a": fa,
                "a_lines": [toks[fa][ia][1], toks[fa][ie + window - 1][1]],
                "b": fb,
                "b_lines": [toks[fb][ib][1], toks[fb][ib + span + window - 1][1]],
            }
        )
    out_regions.sort(key=lambda r: (-r["tokens"], r["a"], r["b"]))
    lines = {name: len({ln for _, ln in seq}) for name, seq in toks.items()}
    files = [
        {"file": name, "duplicated_lines": len(v), "lines": lines[name]}
        for name, v in dup.items()
    ]
    files.sort(key=lambda f: (-f["duplicated_lines"], f["file"]))
    total = sum(lines.values())
    duplicated = sum(len(v) for v in dup.values())
    return {
        "files_scanned": len(toks),
        "lines": total,
        "duplicated_lines": duplicated,
        "density": (100.0 * duplicated / total) if total else 0.0,
        "files": files,
        "regions": out_regions,
    }


def report(lang, result, top, out):
    print(
        f"== {lang}: {result['files_scanned']} files, {result['lines']} code lines, "
        f"{result['duplicated_lines']} duplicated ({result['density']:.1f}%)",
        file=out,
    )
    print("-- files with the most duplicated lines", file=out)
    for f in result["files"][:top]:
        pct = 100.0 * f["duplicated_lines"] / f["lines"] if f["lines"] else 0.0
        print(f"{f['duplicated_lines']:7d} {pct:5.1f}%  {f['file']}", file=out)
    print("-- largest regions shared by two files", file=out)
    for r in result["regions"][:top]:
        print(
            f"{r['tokens']:7d} tok  {r['a']}:{r['a_lines'][0]}-{r['a_lines'][1]}"
            f"  ==  {r['b']}:{r['b_lines'][0]}-{r['b_lines'][1]}",
            file=out,
        )
    print(file=out)


LANGS = {
    "cpp": (CPP_ROOTS, (".cpp", ".hpp", ".h", ".c"), cpp_tokens),
    "python": (PY_ROOTS, (".py",), py_tokens),
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--lang", choices=["cpp", "python", "all"], default="all")
    parser.add_argument(
        "--root", default=REPO, help="repository root (default: this checkout)"
    )
    parser.add_argument(
        "--paths",
        nargs="+",
        help="scan these files or directories instead of the default roots",
    )
    parser.add_argument(
        "--tokens", type=int, default=100, help="window size (default 100)"
    )
    parser.add_argument("--max-occurrences", type=int, default=40)
    parser.add_argument(
        "--top", type=int, default=25, help="rows per list (default 25)"
    )
    parser.add_argument(
        "--json", metavar="PATH", help="also write the full result as JSON"
    )
    args = parser.parse_args(argv)
    langs = ["cpp", "python"] if args.lang == "all" else [args.lang]
    results = {}
    for lang in langs:
        roots, suffixes, tokenizer = LANGS[lang]
        roots = args.paths if args.paths else roots
        paths = []
        for rel in roots:
            full = os.path.join(args.root, rel)
            if os.path.isfile(full):
                if full.endswith(suffixes):
                    paths.append(full)
            else:
                paths.extend(collect(args.root, [rel], suffixes))
        results[lang] = analyse(
            paths, tokenizer, args.tokens, args.max_occurrences, args.root
        )
        report(lang, results[lang], args.top, sys.stdout)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
