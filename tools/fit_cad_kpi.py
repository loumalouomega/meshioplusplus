#!/usr/bin/env python3
"""fit_cad's KPI table (doc/fit_cad_freeform_plan.md, step 1).

Runs ``meshioplusplus_fit_cad_kpi`` (built with the tests:
``build/configure.sh --tests --build``) over the golden, freeform and
kinematic fixtures in every freeform mode and with both backends where the
build has OpenCASCADE, writes the rows as CSV, and regenerates the KPI table
between the ``fit_cad_kpi`` markers of doc/benchmarks.md:

    tools/fit_cad_kpi.py --binary build/cpp-release/meshioplusplus_fit_cad_kpi \\
        --csv benchmark/fit_cad_kpi.csv --doc doc/benchmarks.md
    tools/fit_cad_kpi.py --compare benchmark/fit_cad_kpi.csv --csv new.csv

``--compare OLD`` prints, per fixture, backend and mode, the change of the
faceted area share, of the validity and closedness and of the STEP round trip
against an earlier CSV, and exits 1 when a run lost validity or closedness or
left more of its area as facets: the gate a step of the plan reports before
the next one. ``--option key=value`` and ``--tag`` pass through to the binary,
so two runs can differ in one option; ``--hook PATH`` adds the
topology-optimised hook of notebook 39 (``example/data/hook.mdpa``; minutes to
hours with OpenCASCADE). Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_BINARY = ROOT / "build" / "cpp-release" / "meshioplusplus_fit_cad_kpi"
BEGIN = "<!-- fit_cad_kpi:begin -->"
END = "<!-- fit_cad_kpi:end -->"
MODES = ("legacy", "adaptive", "network", "auto")


def _read(path: pathlib.Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _key(row: dict[str, str]) -> tuple[str, str, str]:
    return (row["fixture"], row["backend"], row["mode"])


def _f(row: dict[str, str], name: str) -> float:
    try:
        return float(row[name])
    except (KeyError, ValueError):
        return math.nan


def _fmt(value: float, digits: int = 3) -> str:
    if math.isnan(value):
        return "–"
    if value == 0.0:
        return "0"
    return f"{value:.{digits}g}"


def run(args: argparse.Namespace) -> list[dict[str, str]]:
    cmd = [str(args.binary)]
    for flag in ("set", "fixtures", "modes", "backends", "hook", "tag"):
        value = getattr(args, flag)
        if value:
            cmd += [f"--{flag}", value]
    if args.variants:
        cmd.append("--variants")
    for opt in args.option or []:
        cmd += ["--option", opt]
    cmd += ["--csv", str(args.csv)]
    print("+", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)
    return _read(args.csv)


def table(rows: list[dict[str, str]]) -> str:
    """The doc table: the faceted area share per fixture, backend and mode,
    then each fixture's validation measures in its best mode."""
    backends = sorted({r["backend"] for r in rows})
    fixtures: list[str] = []
    for r in rows:
        if r["fixture"] not in fixtures:
            fixtures.append(r["fixture"])
    by = {_key(r): r for r in rows}
    out = []
    head = ["Fixture"] + [f"{b} {m}" for b in backends for m in MODES]
    out.append("| " + " | ".join(head) + " |")
    out.append("| " + " | ".join(["---"] * len(head)) + " |")
    for fx in fixtures:
        cells = [fx]
        for b in backends:
            for m in MODES:
                r = by.get((fx, b, m))
                if r is None:
                    cells.append("")
                    continue
                cell = _fmt(_f(r, "facet_area_fraction"))
                if r["valid"] != "1" or r["closed"] != "1":
                    cell += " ✗"
                cells.append(cell)
        out.append("| " + " | ".join(cells) + " |")
    out.append("")
    head = [
        "Fixture",
        "Run",
        "Faces (B-spline)",
        "Control points",
        "Mesh→CAD max / P95",
        "CAD→mesh max / P95",
        "Normal max / P95 (°)",
        "G1 max (°)",
        "STEP read back",
        "Time (s)",
    ]
    out.append("| " + " | ".join(head) + " |")
    out.append("| " + " | ".join(["---"] * len(head)) + " |")
    for fx in fixtures:
        # The run that leaves the least area as facets, then the fastest.
        runs = [r for r in rows if r["fixture"] == fx]
        best = min(runs, key=lambda r: (_f(r, "facet_area_fraction"), _f(r, "seconds")))
        rt = {"1": "valid", "0": "NOT valid", "-1": "–"}.get(
            best["step_round_trip"], "?"
        )
        if best["step_round_trip"] == "1":
            rt += f" ({_fmt(_f(best, 'step_volume_ratio'), 6)})"
        tol = _f(best, "tolerance")
        out.append(
            "| "
            + " | ".join(
                [
                    fx,
                    f"{best['backend']} {best['mode']}",
                    f"{best['num_faces']} ({best['num_bspline']})",
                    best["total_ctrl"],
                    f"{_fmt(_f(best, 'mesh_to_cad_max') / tol)} / "
                    f"{_fmt(_f(best, 'mesh_to_cad_p95') / tol)} tol",
                    f"{_fmt(_f(best, 'cad_to_mesh_max') / tol)} / "
                    f"{_fmt(_f(best, 'cad_to_mesh_p95') / tol)} tol",
                    f"{_fmt(_f(best, 'normal_max_deg'))} / {_fmt(_f(best, 'normal_p95_deg'))}",
                    _fmt(_f(best, "max_g1_deg")),
                    rt,
                    _fmt(_f(best, "seconds")),
                ]
            )
            + " |"
        )
    return "\n".join(out) + "\n"


def write_doc(doc: pathlib.Path, text: str) -> None:
    content = doc.read_text()
    if BEGIN not in content or END not in content:
        raise SystemExit(f"{doc}: no {BEGIN} ... {END} markers")
    head, rest = content.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    doc.write_text(head + BEGIN + "\n\n" + text + "\n" + END + tail)


def compare(old_rows: list[dict[str, str]], new_rows: list[dict[str, str]]) -> int:
    old = {_key(r): r for r in old_rows}
    worse = 0
    print(
        f"{'fixture':24s} {'backend':7s} {'mode':9s} {'facets old':>10s} {'new':>10s}  notes"
    )
    for r in new_rows:
        o = old.get(_key(r))
        if o is None:
            print(
                f"{r['fixture']:24s} {r['backend']:7s} {r['mode']:9s} {'':>10s} "
                f"{_fmt(_f(r, 'facet_area_fraction')):>10s}  new"
            )
            continue
        notes = []
        fo, fn = _f(o, "facet_area_fraction"), _f(r, "facet_area_fraction")
        if fn > fo + 1e-9:
            notes.append("MORE FACETS")
            worse += 1
        for flag in ("valid", "closed"):
            if o[flag] == "1" and r[flag] != "1":
                notes.append(f"LOST {flag.upper()}")
                worse += 1
        if o["step_round_trip"] == "1" and r["step_round_trip"] == "0":
            notes.append("STEP READ-BACK FAILS")
            worse += 1
        if fn < fo - 1e-9:
            notes.append("fewer facets")
        print(
            f"{r['fixture']:24s} {r['backend']:7s} {r['mode']:9s} {_fmt(fo):>10s} "
            f"{_fmt(fn):>10s}  {', '.join(notes)}"
        )
    return 1 if worse else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--binary", type=pathlib.Path, default=DEFAULT_BINARY)
    ap.add_argument(
        "--csv", type=pathlib.Path, default=ROOT / "benchmark" / "fit_cad_kpi.csv"
    )
    ap.add_argument(
        "--doc",
        type=pathlib.Path,
        default=None,
        help="regenerate the table between the markers of this page",
    )
    ap.add_argument(
        "--compare",
        type=pathlib.Path,
        default=None,
        help="an earlier CSV to print the deltas against",
    )
    ap.add_argument("--no-run", action="store_true", help="use --csv as it is")
    ap.add_argument("--set")
    ap.add_argument("--fixtures")
    ap.add_argument("--modes")
    ap.add_argument("--backends")
    ap.add_argument("--hook")
    ap.add_argument("--tag")
    ap.add_argument("--variants", action="store_true")
    ap.add_argument("--option", action="append")
    args = ap.parse_args()
    rows = _read(args.csv) if args.no_run else run(args)
    if args.doc:
        write_doc(args.doc, table(rows))
    if args.compare:
        return compare(_read(args.compare), rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
