#!/usr/bin/env python3
"""Heap allocation counts and peak live bytes per operation row (roadmap 5.28).

Runs every tier S row of ``meshioplusplus_bench_ops`` built SEQ with
``-DMESHIOPLUSPLUS_BENCH_ALLOCS=ON``, each row in its own process on one
thread. The binary's ``--allocs`` mode brackets one post-warmup run with an
opt-in ``operator new``/``delete`` counter (linked into the benchmark binary
only, never the library) that records total allocations, total bytes allocated
and peak live bytes of that run alone.

    tools/bench_allocs.py --record benchmark/allocs_S.csv
    tools/bench_allocs.py --check benchmark/allocs_S.csv

``--check`` flags any row whose ``peak_live_bytes`` rises more than
``--peak-tolerance`` (default 10%) or whose ``allocs`` rises more than
``--alloc-tolerance`` (default 20%) above the committed baseline, or any
missing baseline row. Pass ``--advisory`` to report without failing the
process (used by the pull-request check in ``ci.yml``).

Without ``--binary``, a SEQ Release tree with allocation counting enabled is
configured and built in ``build/bench-allocs`` (CMake and Ninja). Standard
library only.
"""

import argparse
import concurrent.futures
import csv
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def build(tree):
    """Configure and build a SEQ Release ``meshioplusplus_bench_ops_allocs``."""
    subprocess.run(
        [
            "cmake",
            "-S",
            str(ROOT),
            "-B",
            str(tree),
            "-G",
            "Ninja",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DMESHIOPLUSPLUS_BUILD_PYTHON=OFF",
            "-DMESHIOPLUSPLUS_BUILD_BENCHMARKS=ON",
            "-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ",
            "-DMESHIOPLUSPLUS_WITH_HDF5=OFF",
            "-DMESHIOPLUSPLUS_WITH_NETCDF=OFF",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        [
            "cmake",
            "--build",
            str(tree),
            "--target",
            "meshioplusplus_bench_ops_allocs",
        ],
        check=True,
    )
    return tree / "meshioplusplus_bench_ops_allocs"


def toolchain(binary):
    """The compiler (from the tree's CMakeCache.txt)."""
    compiler = os.environ.get("CXX", "c++")
    cache = Path(binary).parent / "CMakeCache.txt"
    if cache.is_file():
        for line in cache.read_text(errors="replace").splitlines():
            if line.startswith("CMAKE_CXX_COMPILER:"):
                compiler = line.split("=", 1)[1]
                break

    def first_line(cmd):
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except (OSError, subprocess.CalledProcessError):
            return "unknown"
        return (out.stdout or out.stderr).splitlines()[0].strip()

    return {"compiler": first_line([compiler, "--version"])}


def _env():
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["MESHIOPLUSPLUS_NUM_THREADS"] = "1"
    return env


def list_rows(binary, tier, ops):
    cmd = [str(binary), "--tier", tier, "--list"]
    if ops:
        cmd += ["--ops", ",".join(ops)]
    out = subprocess.run(
        cmd, capture_output=True, text=True, check=True, env=_env()
    ).stdout
    return [r for r in out.split() if r]


def count(binary, tier, op):
    """(cells, allocs, total_bytes, peak_live_bytes) of one row, or None."""
    proc = subprocess.run(
        [
            str(binary),
            "--tier",
            tier,
            "--ops",
            op,
            "--allocs",
        ],
        capture_output=True,
        text=True,
        env=_env(),
    )
    if proc.returncode != 0:
        raise RuntimeError(f"{op}: exit {proc.returncode}\n{proc.stderr[-2000:]}")
    rows = [r for r in csv.DictReader(proc.stdout.splitlines()) if r["op"] == op]
    if not rows:
        return None
    row = rows[0]
    return (
        int(row["cells"]),
        int(row["allocs"]),
        int(row["total_bytes"]),
        int(row["peak_live_bytes"]),
    )


def measure(binary, tier, ops, jobs):
    rows = list_rows(binary, tier, ops)
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(count, binary, tier, op): op for op in rows}
        for future in concurrent.futures.as_completed(futures):
            op = futures[future]
            got = future.result()
            if got is None:
                print(f"skipped by the binary: {op}", file=sys.stderr)
                continue
            results[op] = got
            print(
                f"{op}: allocs={got[1]:,} total={got[2]:,}B peak={got[3]:,}B",
                file=sys.stderr,
            )
    return {op: results[op] for op in rows if op in results}


def write_csv(path, tier, tools, results):
    with open(path, "w", newline="") as fh:
        fh.write(f"# tier: {tier}\n")
        fh.write(f"# compiler: {tools['compiler']}\n")
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["op", "cells", "allocs", "total_bytes", "peak_live_bytes"])
        for op, (cells, allocs, total_bytes, peak_live) in results.items():
            writer.writerow([op, cells, allocs, total_bytes, peak_live])


def read_csv(path):
    meta, lines = {}, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            key, _, value = line[1:].partition(":")
            meta[key.strip()] = value.strip()
        elif line.strip():
            lines.append(line)
    rows = {
        r["op"]: (
            int(r["cells"]),
            int(r["allocs"]),
            int(r["total_bytes"]),
            int(r["peak_live_bytes"]),
        )
        for r in csv.DictReader(lines)
    }
    return meta, rows


def _ratio(now_val, base_val):
    if base_val == 0:
        return 1.0 if now_val == 0 else float("inf")
    return now_val / base_val


def compare(baseline, results, peak_tolerance, alloc_tolerance):
    """Compare measured rows against baseline; return (report, failed)."""
    report, failed = [], False
    for op in sorted(set(baseline) | set(results)):
        base = baseline.get(op)
        now = results.get(op)
        if now is None:
            report.append((op, base, None, None, None, "MISSING"))
            failed = True
            continue
        if base is None:
            report.append((op, None, now, None, None, "new row (not gated)"))
            continue
        if base[0] != now[0]:
            report.append((op, base, now, None, None, "CELLS CHANGED"))
            failed = True
            continue
        r_alloc = _ratio(now[1], base[1])
        r_peak = _ratio(now[3], base[3])
        flags = []
        if r_peak > 1.0 + peak_tolerance:
            flags.append(f"PEAK +{r_peak - 1.0:.1%}")
            failed = True
        if r_alloc > 1.0 + alloc_tolerance:
            flags.append(f"ALLOCS +{r_alloc - 1.0:.1%}")
            failed = True
        if not flags:
            if r_peak < 1.0 - peak_tolerance or r_alloc < 1.0 - alloc_tolerance:
                verdict = "smaller: re-record the baseline"
            else:
                verdict = ""
        else:
            verdict = ", ".join(flags)
        report.append((op, base, now, r_alloc, r_peak, verdict))
    return report, failed


def markdown(report, peak_tolerance, alloc_tolerance, note):
    out = ["### Heap allocations and peak live memory, tier S, one thread\n"]
    if note:
        out.append(note + "\n")
    out.append(
        f"Tolerances: peak live bytes +{peak_tolerance:.0%}, "
        f"allocation count +{alloc_tolerance:.0%}. "
        "Rows within both tolerances are not listed.\n"
    )
    out.append(
        "| row | allocs (base -> now) | ratio | peak live bytes (base -> now) | ratio | |"
    )
    out.append("| --- | ---: | ---: | ---: | ---: | --- |")
    shown = 0
    for op, base, now, r_alloc, r_peak, verdict in report:
        if not verdict:
            continue
        shown += 1
        a_col = (
            f"{base[1]:,} -> {now[1]:,}"
            if base is not None and now is not None
            else (f"{base[1]:,} ->" if base is not None else f"-> {now[1]:,}")
        )
        p_col = (
            f"{base[3]:,} -> {now[3]:,}"
            if base is not None and now is not None
            else (f"{base[3]:,} ->" if base is not None else f"-> {now[3]:,}")
        )
        ra = "" if r_alloc is None else f"{r_alloc:.3f}"
        rp = "" if r_peak is None else f"{r_peak:.3f}"
        out.append(f"| `{op}` | {a_col} | {ra} | {p_col} | {rp} | {verdict} |")
    if not shown:
        out.append("| every row within tolerances | | | | | |")
    return "\n".join(out) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", metavar="CSV", help="write counts as a baseline")
    mode.add_argument("--check", metavar="CSV", help="compare against this baseline")
    parser.add_argument("--binary", help="a SEQ meshioplusplus_bench_ops_allocs to use")
    parser.add_argument(
        "--build-dir",
        default=str(ROOT / "build" / "bench-instructions"),
        help="where to build one without --binary",
    )
    parser.add_argument("--tier", default="S")
    parser.add_argument("--ops", help="comma-separated rows (default: every row)")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--peak-tolerance", type=float, default=0.10)
    parser.add_argument("--alloc-tolerance", type=float, default=0.20)
    parser.add_argument(
        "--advisory",
        action="store_true",
        help="report regressions without returning a non-zero exit status",
    )
    parser.add_argument("--out", help="with --check: also write this run's counts")
    parser.add_argument("--summary", help="with --check: append a Markdown report")
    args = parser.parse_args(argv)

    binary = Path(args.binary) if args.binary else build(Path(args.build_dir))
    ops = [o for o in (args.ops or "").split(",") if o]
    tools = toolchain(binary)
    results = measure(binary, args.tier, ops, args.jobs)

    if args.record:
        write_csv(args.record, args.tier, tools, results)
        print(f"wrote {args.record}: {len(results)} rows", file=sys.stderr)
        return 0

    meta, baseline = read_csv(args.check)
    if ops:
        baseline = {op: v for op, v in baseline.items() if op in ops}
    if args.out:
        write_csv(args.out, args.tier, tools, results)
    report, failed = compare(
        baseline, results, args.peak_tolerance, args.alloc_tolerance
    )
    note = ""
    same_tools = meta.get("compiler") == tools["compiler"]
    if not same_tools:
        note = (
            f"The baseline was recorded with `{meta.get('compiler')}`, this run "
            f"with `{tools['compiler']}`: the counts are not directly comparable."
        )
    text = markdown(report, args.peak_tolerance, args.alloc_tolerance, note)
    print(text)
    if args.summary:
        with open(args.summary, "a") as fh:
            fh.write(text)
    if args.advisory:
        return 0
    return 1 if failed and same_tools else 0


if __name__ == "__main__":
    sys.exit(main())
