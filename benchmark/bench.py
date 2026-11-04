"""Benchmark harness: meshio++ read/write timings, optionally against legacy meshio.

Two modes share the timing code:

* the curated comparison the notebook plots (``FORMATS``: a handful of formats,
  legacy pure-Python ``meshio`` against meshio++);
* a sweep over **every** readable and writable format (``all_formats()``,
  roadmap §3), each fed the synthetic mesh its conformance declaration says it
  keeps (``tests/python/conformance_spec.py``): tetrahedra where it takes them,
  the tetrahedra's surface for a surface format, the points for a cloud.

Legacy meshio is optional. It is imported from ``MESHIO_LEGACY_SRC`` (a
source checkout's ``src`` directory; it is pure Python) or as an installed
``meshio``; without either, only the meshio++ columns are filled.

Usage (see also 01_benchmark.ipynb)::

    python benchmark/bench.py --sizes S,M --formats vtu,gmsh --out results_all.csv
    python benchmark/bench.py --sizes L --formats vtu,stl_binary --memory

or, from Python::

    from bench import FORMATS, run
    records = run(points, cells, FORMATS, repeats=3)
"""

from __future__ import annotations

import argparse
import os
import pathlib
import statistics
import sys
import tempfile
import time

import numpy as np

import meshioplusplus as pp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "tests" / "python"))


def _import_legacy():
    src = os.environ.get("MESHIO_LEGACY_SRC")
    if src and src not in sys.path:
        sys.path.insert(0, src)
    try:
        import meshio
    except ImportError:
        return None, None
    version = getattr(meshio, "__version__", "?")
    if src:
        try:
            import tomllib

            with open(os.path.join(src, "..", "pyproject.toml"), "rb") as fh:
                version = tomllib.load(fh)["project"]["version"]
        except (OSError, KeyError, ImportError):
            pass
    return meshio, version


legacy, LEGACY_VERSION = _import_legacy()
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

# Grid points per edge of the synthetic tetrahedral cube, per size tier
# (6 (n-1)^3 tetrahedra: about 20k, 250k and 1M).
SIZES = {"S": 16, "M": 36, "L": 56}

# These formats reject the mixed conformance mesh, but accept the harness's
# single-type input. Keep this explicit rather than treating unknown ("?")
# conformance cells as proven round trips.
SINGLE_TYPE_INPUTS = {"mfm": "volume"}
POINT_FIELD_INPUTS = {"dex": "benchmark:field"}

# Rows for an encoding a format's default write does not reach: name ->
# (format, write keyword arguments, input kind). STL writes ASCII by default.
VARIANTS = {"stl_binary": ("stl", {"binary": True}, "surface")}


def all_formats():
    """One spec per format meshio++ both writes and reads back, by input kind.

    Returns ``[(format, kind)]`` with kind ``"volume"``, ``"surface"`` or
    ``"points"`` -- the largest thing the format's conformance declaration
    says survives a round trip, plus explicitly tested single-type exceptions
    in ``SINGLE_TYPE_INPUTS`` and required nodal-field inputs in
    ``POINT_FIELD_INPUTS``. Other failing, write-only or entirely lossy
    formats are left out.
    """
    import conformance_spec as cs

    out = []
    for fmt, spec in sorted(cs.SPEC.items()):
        if fmt in POINT_FIELD_INPUTS:
            out.append((fmt, "points"))
            continue
        if fmt in SINGLE_TYPE_INPUTS:
            out.append((fmt, SINGLE_TYPE_INPUTS[fmt]))
            continue
        if "points" not in spec:
            continue
        kept = {t for t, o in spec["cells"].items() if o in ("exact", "reordered")}
        if "tetra" in kept:
            out.append((fmt, "volume"))
        elif "triangle" in spec["cells"] and spec.get("input") != "2d":
            # Accepted on its own (STL loses triangles only next to volumes).
            out.append((fmt, "surface"))
        elif spec["points"] in ("exact", "approx"):
            out.append((fmt, "points"))
    out += [(name, kind) for name, (_, _, kind) in VARIANTS.items()]
    return sorted(out)


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
    path = pathlib.Path(path)
    if path.is_dir():
        return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())
    # A format that writes siblings (xdmf + h5, tetgen .node + .ele, ...).
    return sum(
        f.stat().st_size for f in path.parent.glob(path.stem + ".*") if f.is_file()
    )


def bench_one(points, cells, spec, tmpdir, repeats):
    """Benchmark one curated format; returns a record dict (or None if unsupported)."""
    label, fn, wkw, rfmt, cpp = spec
    mesh_pp = pp.Mesh(points, cells)
    pp_path = os.path.join(tmpdir, "pp_" + fn)
    try:
        mesh_pp.write(pp_path, **wkw)
        pp.read(pp_path, file_format=rfmt)
    except Exception as exc:  # a format this build cannot do: reported, skipped
        print(f"  skip {label}: {type(exc).__name__}: {exc}")
        return None
    rec = {
        "format": label,
        "cpp": cpp,
        "cells": sum(len(c[1]) for c in cells),
        "bytes": _total_size(pp_path),
        "write_pp": _median(lambda: mesh_pp.write(pp_path, **wkw), repeats),
        "read_pp": _median(lambda: pp.read(pp_path, file_format=rfmt), repeats),
        "write_legacy": float("nan"),
        "read_legacy": float("nan"),
    }
    if legacy is not None:
        mesh_lg = legacy.Mesh(points, cells)
        lg_path = os.path.join(tmpdir, "lg_" + fn)
        try:
            mesh_lg.write(lg_path, **wkw)
            legacy.read(pp_path, file_format=rfmt)
        except Exception as exc:  # legacy lacks the format: meshio++ columns only
            print(f"  legacy skips {label}: {type(exc).__name__}")
        else:
            rec["write_legacy"] = _median(
                lambda: mesh_lg.write(lg_path, **wkw), repeats
            )
            # Both read the *same* reference file (written by meshio++).
            rec["read_legacy"] = _median(
                lambda: legacy.read(pp_path, file_format=rfmt), repeats
            )
    rec["write_speedup"] = rec["write_legacy"] / rec["write_pp"]
    rec["read_speedup"] = rec["read_legacy"] / rec["read_pp"]
    return rec


def run(points, cells, formats=FORMATS, repeats=3):
    """Benchmark every curated format; returns a list of record dicts."""
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


def _inputs(size):
    from inputs import synthetic_tet_grid

    points, cells = synthetic_tet_grid(SIZES[size])
    volume = pp.Mesh(points, cells)
    surface = pp.extract_surface(volume)
    return {
        "volume": volume,
        "surface": pp.Mesh(surface.points, [(b.type, b.data) for b in surface.cells]),
        "points": pp.Mesh(points, [("vertex", np.arange(len(points))[:, None])]),
    }


def _input_for_format(mesh, fmt):
    """Supply required fields without mutating shared geometry-only inputs."""
    if fmt not in POINT_FIELD_INPUTS:
        return mesh
    return pp.Mesh(
        mesh.points,
        [(block.type, block.data) for block in mesh.cells],
        point_data={
            POINT_FIELD_INPUTS[fmt]: np.arange(len(mesh.points), dtype=np.float64)
        },
    )


def _status_mb(field):
    """``VmRSS`` or ``VmHWM`` from ``/proc/self/status``, in MB (Linux)."""
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith(field + ":"):
                return int(line.split()[1]) / 1024
    raise OSError(field)


def _reset_peak():
    """Start a fresh high-water mark: the RSS now, with ``VmHWM`` reset to it
    where the kernel allows (Linux ``clear_refs`` 5); elsewhere the process's
    peak so far (``ru_maxrss``, KB on Linux and bytes on macOS), which the
    measurement then only exceeds."""
    try:
        with open("/proc/self/clear_refs", "w") as fh:
            fh.write("5")
        return _status_mb("VmRSS")
    except OSError:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / (1 << 20) if sys.platform == "darwin" else peak / 1024


def _peak_mb():
    try:
        return _status_mb("VmHWM")
    except OSError:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / (1 << 20) if sys.platform == "darwin" else peak / 1024


def _save_input(mesh, path):
    """The input geometry as an ``.npz`` a fresh process loads without the
    transient peak of building it (``extract_surface`` and the like)."""
    arrays = {"points": mesh.points}
    for k, block in enumerate(mesh.cells):
        arrays[f"cells{k}_{block.type}"] = block.data
    for name, values in mesh.point_data.items():
        arrays[f"point_data:{name}"] = values
    np.savez(path, **arrays)


def _load_input(path):
    with np.load(path) as z:
        cells = [
            (key.split("_", 1)[1], z[key]) for key in z.files if key.startswith("cells")
        ]
        point_data = {
            key.split(":", 1)[1]: z[key]
            for key in z.files
            if key.startswith("point_data:")
        }
        return pp.Mesh(z["points"], cells, point_data=point_data)


def _measure(argv):
    """The child of ``--memory``: one write or one read, then the peak RSS it
    added over what the process held before it, printed in MB."""
    op, fmt, path, input_path, write_kwargs = argv[0], argv[1], argv[2], argv[3], {}
    import json

    write_kwargs = json.loads(argv[4])
    mesh = _load_input(input_path) if op == "write" else None
    # Load the format and the core first: the measurement is the operation's.
    from meshioplusplus import _core  # noqa: F401
    from meshioplusplus._helpers import _ensure_format

    _ensure_format(fmt)
    before = _reset_peak()
    if op == "write":
        pp.write(path, mesh, file_format=fmt, **write_kwargs)
    else:
        pp.read(path, file_format=fmt)
    print(f"{_peak_mb() - before:.1f}")
    return 0


def _peak_in_child(op, fmt, path, input_path, write_kwargs):
    """Peak RSS (MB) one write or read adds, in a fresh interpreter, so one
    row's high-water mark cannot hide the next one's."""
    import json
    import subprocess

    out = subprocess.run(
        [sys.executable, __file__, "--_measure", op, fmt, str(path), str(input_path)]
        + [json.dumps(write_kwargs)],
        check=True,
        capture_output=True,
        text=True,
    )
    return float(out.stdout.strip().splitlines()[-1])


def run_all(sizes=("S",), formats=None, repeats=3, memory=False):
    """Time a write and a read of every registry format; returns records.

    With ``memory``, each row also carries the peak resident memory its write
    and its read add (``write_peak_mb``, ``read_peak_mb``), measured in a
    fresh process per operation (roadmap 5.1.3)."""
    import conformance_spec as cs

    wanted = set(formats) if formats else None
    records = []
    for size in sizes:
        meshes = _inputs(size)
        for name, kind in all_formats():
            if wanted and name not in wanted:
                continue
            fmt, write_kwargs, _ = VARIANTS.get(name, (name, {}, kind))
            mesh = _input_for_format(meshes[kind], fmt)
            with tempfile.TemporaryDirectory() as tmp:
                path = cs._target(pathlib.Path(tmp), fmt)
                if path.suffix == "":
                    path.mkdir(parents=True, exist_ok=True)
                read_as = cs._READ_AS.get(fmt, fmt)
                try:
                    pp.write(path, mesh, file_format=fmt, **write_kwargs)
                    pp.read(path, file_format=read_as)
                except Exception as exc:  # reported, not fatal: the sweep goes on
                    print(f"  skip {name} ({size}): {type(exc).__name__}: {exc}")
                    continue
                rec = {
                    "size": size,
                    "format": name,
                    "input": kind,
                    "cells": sum(len(b.data) for b in mesh.cells),
                    "bytes": _total_size(path),
                    "write_s": _median(
                        lambda: pp.write(path, mesh, file_format=fmt, **write_kwargs),
                        repeats,
                    ),
                    "read_s": _median(
                        lambda: pp.read(path, file_format=read_as), repeats
                    ),
                }
                if memory:
                    saved = pathlib.Path(tmp) / "input.npz"
                    _save_input(mesh, saved)
                    rec["write_peak_mb"] = _peak_in_child(
                        "write", fmt, path, saved, write_kwargs
                    )
                    rec["read_peak_mb"] = _peak_in_child(
                        "read", read_as, path, saved, {}
                    )
            records.append(rec)
            peaks = (
                f"  peak +{rec['write_peak_mb']:7.1f} / +{rec['read_peak_mb']:7.1f} MB"
                if memory
                else ""
            )
            print(
                f"  {size} {name:14s} {kind:8s} write {rec['write_s'] * 1e3:8.1f} ms  "
                f"read {rec['read_s'] * 1e3:8.1f} ms{peaks}"
            )
    return records


def write_csv(records, path):
    import csv

    fields = list(records[0]) if records else []
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in records:
            w.writerow({k: r[k] for k in fields})


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Time meshio++ read/write for every format"
    )
    ap.add_argument("--sizes", default="S", help="comma-separated tiers: S, M, L")
    ap.add_argument("--formats", default="", help="comma-separated subset")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--out", default="results_all.csv")
    ap.add_argument(
        "--memory",
        action="store_true",
        help="also record each write's and read's peak RSS, one subprocess each",
    )
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["--_measure"]:
        return _measure(argv[1:])
    args = ap.parse_args(argv)
    sizes = [s for s in args.sizes.split(",") if s]
    unknown = set(sizes) - set(SIZES)
    if unknown:
        ap.error(f"unknown size(s): {', '.join(sorted(unknown))}")
    records = run_all(
        sizes, [f for f in args.formats.split(",") if f], args.repeats, args.memory
    )
    write_csv(records, args.out)
    print(f"wrote {len(records)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
