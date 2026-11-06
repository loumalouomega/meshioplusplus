#!/usr/bin/env python3
"""Instruction counts of the operation rows, and a gate on them (roadmap 5.1.1).

Timings never fail a build: hosted runners are noisy (doc/benchmark_trends.md).
An instruction count is stable where a wall time is not, so this counts the
instructions of one run of every tier S row of ``meshioplusplus_bench_ops``,
built SEQ and run on one thread, each row in its own process under
``valgrind --tool=cachegrind --instr-at-start=no``: the binary's
``--instructions`` mode brackets the row's run with Cachegrind's start and
stop requests, after the untimed warmup that builds the row's inputs, so the
count is the operation's alone.

    tools/bench_instructions.py --record benchmark/instructions_S.csv
    tools/bench_instructions.py --check benchmark/instructions_S.csv

``--check`` fails when a row counts more than ``--tolerance`` (default 2%)
above the committed baseline, or when a baseline row is missing; a row that
counts less than the tolerance below its baseline, or a new row, is reported
but passes, and the baseline should then be re-recorded. A count depends on
the compiler, its flags and the C++ library, so the baseline names the
toolchain it was recorded with: against another, the comparison is printed
and nothing fails. A change that moves a count on purpose re-records the
baseline in the same pull request and says why (doc/benchmarks.md,
"Instruction-count gate").

Without ``--binary``, a SEQ Release tree is configured and built in
``build/bench-instructions`` (CMake and Ninja, as tools/bench_ops.sh does).
Standard library only.
"""

import argparse
import concurrent.futures
import csv
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Rows left out of the gate, each with its reason. A row costs about thirty
# times its native time under Cachegrind; these would hold the job for minutes
# each and are covered by the weekly timings instead.
SKIP = {}

_IREFS = re.compile(r"==\d+==\s+I\s+refs:\s+([\d,]+)")


def build(tree):
    """Configure and build a SEQ Release ``meshioplusplus_bench_ops``."""
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
        ["cmake", "--build", str(tree), "--target", "meshioplusplus_bench_ops"],
        check=True,
    )
    return tree / "meshioplusplus_bench_ops"


def toolchain(binary):
    """The compiler (from the tree's CMakeCache.txt) and valgrind versions."""
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

    return {
        "compiler": first_line([compiler, "--version"]),
        "valgrind": first_line(["valgrind", "--version"]),
    }


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
    """(cells, instructions) of one row, or None when the binary skipped it."""
    proc = subprocess.run(
        [
            "valgrind",
            "--tool=cachegrind",
            "--cache-sim=no",
            "--instr-at-start=no",
            "--cachegrind-out-file=/dev/null",
            str(binary),
            "--tier",
            tier,
            "--ops",
            op,
            "--instructions",
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
    match = _IREFS.search(proc.stderr)
    if not match:
        raise RuntimeError(f"{op}: no 'I refs' in valgrind's output")
    return int(rows[0]["cells"]), int(match.group(1).replace(",", ""))


def measure(binary, tier, ops, jobs):
    rows = [r for r in list_rows(binary, tier, ops) if r not in SKIP]
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
            print(f"{op}: {got[1]:,}", file=sys.stderr)
    # The binary's row order, which is stable and groups related rows.
    return {op: results[op] for op in rows if op in results}


def write_csv(path, tier, tools, results):
    with open(path, "w", newline="") as fh:
        fh.write(f"# tier: {tier}\n")
        fh.write(f"# compiler: {tools['compiler']}\n")
        fh.write(f"# valgrind: {tools['valgrind']}\n")
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["op", "cells", "instructions"])
        for op, (cells, n) in results.items():
            writer.writerow([op, cells, n])


def read_csv(path):
    meta, lines = {}, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            key, _, value = line[1:].partition(":")
            meta[key.strip()] = value.strip()
        elif line.strip():
            lines.append(line)
    rows = {
        r["op"]: (int(r["cells"]), int(r["instructions"]))
        for r in csv.DictReader(lines)
    }
    return meta, rows


def compare(baseline, results, tolerance):
    """The report rows and whether the gate fails: (op, base, now, ratio, verdict)."""
    report, failed = [], False
    for op in sorted(set(baseline) | set(results)):
        base = baseline.get(op)
        now = results.get(op)
        if now is None:
            report.append((op, base[1], None, None, "MISSING"))
            failed = True
            continue
        if base is None:
            report.append((op, None, now[1], None, "new row (not gated)"))
            continue
        if base[0] != now[0]:
            report.append((op, base[1], now[1], None, "CELLS CHANGED"))
            failed = True
            continue
        ratio = now[1] / base[1] if base[1] else float("inf")
        if ratio > 1 + tolerance:
            verdict = "SLOWER"
            failed = True
        elif ratio < 1 - tolerance:
            verdict = "faster: re-record the baseline"
        else:
            verdict = ""
        report.append((op, base[1], now[1], ratio, verdict))
    return report, failed


def markdown(report, tolerance, note):
    out = ["### Instruction counts, tier S, one thread\n"]
    if note:
        out.append(note + "\n")
    out.append(f"Tolerance: {tolerance:.1%}. Rows within it are not listed.\n")
    out.append("| row | baseline | now | ratio | |")
    out.append("| --- | ---: | ---: | ---: | --- |")
    shown = 0
    for op, base, now, ratio, verdict in report:
        if not verdict:
            continue
        shown += 1
        out.append(
            f"| `{op}` | {'' if base is None else f'{base:,}'} | "
            f"{'' if now is None else f'{now:,}'} | "
            f"{'' if ratio is None else f'{ratio:.4f}'} | {verdict} |"
        )
    if not shown:
        out.append("| every row within the tolerance | | | | |")
    return "\n".join(out) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", metavar="CSV", help="write the counts as a baseline")
    mode.add_argument("--check", metavar="CSV", help="compare against this baseline")
    parser.add_argument("--binary", help="a SEQ meshioplusplus_bench_ops to use")
    parser.add_argument(
        "--build-dir",
        default=str(ROOT / "build" / "bench-instructions"),
        help="where to build one without --binary",
    )
    parser.add_argument("--tier", default="S")
    parser.add_argument("--ops", help="comma-separated rows (default: every row)")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--tolerance", type=float, default=0.02)
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
    report, failed = compare(baseline, results, args.tolerance)
    note = ""
    same_tools = meta.get("compiler") == tools["compiler"]
    if not same_tools:
        note = (
            f"The baseline was recorded with `{meta.get('compiler')}`, this run "
            f"with `{tools['compiler']}`: the counts are not comparable, so "
            "nothing fails. Re-record the baseline with this toolchain."
        )
    text = markdown(report, args.tolerance, note)
    print(text)
    if args.summary:
        with open(args.summary, "a") as fh:
            fh.write(text)
    return 1 if failed and same_tools else 0


if __name__ == "__main__":
    sys.exit(main())
