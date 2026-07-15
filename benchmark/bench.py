"""Benchmark harness: legacy pure-Python meshio vs meshio++.

Imports *both* libraries into one process -- ``meshioplusplus`` (the installed
C++-accelerated fork) and the legacy pure-Python ``meshio`` from source, placed
on ``sys.path`` (it is pure Python, no build step). Their ``Mesh``/``read``/
``write`` APIs are identical, so the same geometry is written/read by each and
timed.

Usage (see 01_benchmark.ipynb)::

    from bench import LEGACY_VERSION, FORMATS, run
    records = run(points, cells, FORMATS, repeats=3)
"""

from __future__ import annotations

import os
import statistics
import sys
import tempfile
import time

# --- make the legacy pure-Python meshio importable as ``meshio`` -------------
LEGACY_SRC = "/home/vicente/src/meshio_legacy/src"
if LEGACY_SRC not in sys.path:
    sys.path.insert(0, LEGACY_SRC)

import meshio as legacy  # noqa: E402  (import after sys.path tweak)

import meshioplusplus as pp  # noqa: E402


def _legacy_version() -> str:
    """Read the legacy version from its pyproject (it runs uninstalled)."""
    try:
        import tomllib

        with open(os.path.join(LEGACY_SRC, "..", "pyproject.toml"), "rb") as fh:
            return tomllib.load(fh)["project"]["version"]
    except Exception:
        return legacy.__version__


LEGACY_VERSION = _legacy_version()
PP_VERSION = pp.__version__


# --- format specifications ----------------------------------------------------
# label, filename, write kwargs, read file_format, cpp (does meshio++ use its
# C++ path for this format?)
FORMATS = [
    ("vtu (binary+zlib)", "mesh.vtu", dict(binary=True), None, True),
    ("vtu (ascii)", "mesh_ascii.vtu", dict(binary=False), None, True),
    ("vtk (binary)", "mesh.vtk", dict(binary=True), None, True),
    ("gmsh (binary)", "mesh.msh", dict(binary=True, file_format="gmsh"), "gmsh", True),
    ("xdmf (HDF5)", "mesh.xdmf", {}, None, True),
    ("med (HDF5)", "mesh.med", {}, None, True),
    ("mdpa", "mesh.mdpa", dict(file_format="mdpa"), "mdpa", False),
]


def _median(fn, repeats, warmup=1):
    for _ in range(warmup):
        fn()
    ts = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return statistics.median(ts)


def _total_size(path):
    size = os.path.getsize(path)
    h5 = os.path.splitext(path)[0] + ".h5"
    if path.endswith(".xdmf") and os.path.exists(h5):
        size += os.path.getsize(h5)
    return size


def bench_one(points, cells, spec, tmpdir, repeats):
    """Benchmark one format; returns a record dict (or None if unsupported)."""
    label, fn, wkw, rfmt, cpp = spec
    mesh_pp = pp.Mesh(points, cells)
    mesh_lg = legacy.Mesh(points, cells)

    pp_path = os.path.join(tmpdir, "pp_" + fn)
    lg_path = os.path.join(tmpdir, "lg_" + fn)

    # Smoke-test both directions once; skip the format if either lib can't do it.
    try:
        mesh_pp.write(pp_path, **wkw)
        pp.read(pp_path, file_format=rfmt)
        mesh_lg.write(lg_path, **wkw)
        legacy.read(lg_path, file_format=rfmt)
    except Exception as exc:
        print(f"  skip {label}: {type(exc).__name__}: {exc}")
        return None

    w_pp = _median(lambda: mesh_pp.write(pp_path, **wkw), repeats)
    w_lg = _median(lambda: mesh_lg.write(lg_path, **wkw), repeats)

    # Read timing: both read the *same* reference file (written by meshio++).
    r_pp = _median(lambda: pp.read(pp_path, file_format=rfmt), repeats)
    r_lg = _median(lambda: legacy.read(pp_path, file_format=rfmt), repeats)

    return {
        "format": label,
        "cpp": cpp,
        "bytes": _total_size(pp_path),
        "write_pp": w_pp,
        "write_legacy": w_lg,
        "read_pp": r_pp,
        "read_legacy": r_lg,
        "write_speedup": w_lg / w_pp if w_pp else float("nan"),
        "read_speedup": r_lg / r_pp if r_pp else float("nan"),
    }


def run(points, cells, formats=FORMATS, repeats=3):
    """Benchmark every format; returns a list of record dicts."""
    records = []
    with tempfile.TemporaryDirectory() as tmpdir:
        for spec in formats:
            rec = bench_one(points, cells, spec, tmpdir, repeats)
            if rec is not None:
                records.append(rec)
                print(
                    f"  {rec['format']:20s} "
                    f"write {rec['write_legacy'] * 1e3:7.1f}"
                    f"->{rec['write_pp'] * 1e3:7.1f} ms "
                    f"({rec['write_speedup']:4.1f}x)  "
                    f"read {rec['read_legacy'] * 1e3:7.1f}"
                    f"->{rec['read_pp'] * 1e3:7.1f} ms "
                    f"({rec['read_speedup']:4.1f}x)"
                )
    return records


def write_csv(records, path):
    import csv

    fields = [
        "format",
        "cpp",
        "bytes",
        "write_legacy",
        "write_pp",
        "write_speedup",
        "read_legacy",
        "read_pp",
        "read_speedup",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in records:
            w.writerow({k: r[k] for k in fields})
