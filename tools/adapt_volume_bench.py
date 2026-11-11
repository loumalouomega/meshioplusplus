#!/usr/bin/env python3
"""Benchmark adapt_volume on octree_mesh's tetrahedra (doc/adapt_volume.md).

Each skin is meshed by ``octree_mesh`` at a few ``min_size`` values, then
remeshed by ``adapt_volume`` with sizes matched to the octree's (``hmin`` the
octree's ``min_size``, ``hmax`` its ``max_size``) and a gradation. For each run
the table reports the tetrahedra and points before and after, the dihedral
range, the share of tetrahedra under 10 and 20 degrees, the mean ratio and
scaled-Jacobian minima, the symmetric Hausdorff distance to the skin, the
volume error against the octree mesh, and the seconds.

Baselines, each skipped when absent: MMG's ``adapt`` (the optional
``meshioplusplus-mmg`` add-on, or a build with it) at the same ``hmin``,
``hmax``, ``hgrad`` and ``hausd``, or else the ``mmg3d`` program through Medit
files (``--mmg3d PATH``, or ``mmg3d_O3``/``mmg3d`` on ``PATH``; nothing links
it), and the tetgen and gmsh programs through ``tetrahedralize`` at
``max_size`` = ``hmax``. The classic scans (Happy Buddha,
XYZ RGB dragon) are fetched from common-3d-test-models with ``--fetch``.

Usage::

    python tools/adapt_volume_bench.py [--csv benchmark/adapt_volume.csv]
        [--hgrad 1.3] [--fetch] [--cases sphere,box,cylinder,bunny]
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request

import numpy as np

import meshioplusplus as mp

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLASSICS = (
    "https://raw.githubusercontent.com/alecjacobson/common-3d-test-models/master/data/"
)


def _skin_box():
    return mp.extract_surface(mp.box(dims=(1, 1, 1)))


def _skin_cylinder():
    return mp.extract_surface(mp.cylinder(radius=0.5, height=2.0, nxy=24, nz=8))


def _skin_bunny():
    path = os.path.join(REPO, "example", "data", "Bunny.stl")
    return mp.read(path) if os.path.exists(path) else None


def _skin_classic(name):
    def load():
        path = os.path.join(tempfile.gettempdir(), name)
        try:
            if not os.path.exists(path):
                urllib.request.urlretrieve(CLASSICS + name, path)
            return mp.read(path)
        except Exception as exc:  # pragma: no cover - needs the network
            print(f"{name}: not available ({exc.__class__.__name__}); skipped")
            return None

    return load


# name -> (skin factory, min_size values as fractions of the bounding-box diagonal)
CASES = {
    "sphere": (lambda: mp.icosa_sphere(1.0, 5), (0.02, 0.01)),
    "box": (_skin_box, (0.03, 0.015)),
    "cylinder": (_skin_cylinder, (0.02, 0.01)),
    "bunny": (_skin_bunny, (0.012, 0.008)),
    "buddha": (_skin_classic("happy.obj"), (0.008,)),
    "dragon": (_skin_classic("xyzrgb_dragon.obj"), (0.008,)),
}


def _diagonal(mesh):
    p = np.asarray(mesh.points, dtype=float)
    return float(np.linalg.norm(p.max(axis=0) - p.min(axis=0)))


def _measures(mesh, skin, volume_ref):
    q = mp.attach_quality(mesh)
    dmin = np.concatenate([np.asarray(a) for a in q.cell_data["quality:min_dihedral"]])
    dmax = np.concatenate([np.asarray(a) for a in q.cell_data["quality:max_dihedral"]])
    sj = np.concatenate([np.asarray(a) for a in q.cell_data["quality:scaled_jacobian"]])
    vol = float(
        np.concatenate([np.asarray(a) for a in q.cell_data["quality:volume"]]).sum()
    )
    h = mp.hausdorff_distance(mp.extract_surface(mesh), skin, face_samples=2)
    ntets = int(sum(len(b.data) for b in mesh.cells))
    return {
        "tets": ntets,
        "points": len(mesh.points),
        "min_dihedral": float(dmin.min()),
        "max_dihedral": float(dmax.max()),
        "under_10": float((dmin < 10.0).mean()),
        "under_20": float((dmin < 20.0).mean()),
        "min_scaled_jacobian": float(sj.min()),
        "dihedral_histogram": np.histogram(
            dmin, bins=(0, 5, 10, 20, 30, 40, 50, 60, 90)
        )[0],
        "hausdorff": float(h["distance"]),
        "volume_error": abs(vol - volume_ref) / volume_ref if volume_ref else 0.0,
    }


def _run(rows, case, frac, method, mesh, skin, volume_ref, seconds, extra=""):
    m = _measures(mesh, skin, volume_ref)
    row = {"case": case, "min_size_frac": frac, "method": method, "seconds": seconds}
    row.update({k: v for k, v in m.items() if k != "dihedral_histogram"})
    row["dihedral_histogram"] = " ".join(str(int(x)) for x in m["dihedral_histogram"])
    row["note"] = extra
    rows.append(row)
    print(
        f"{case:9s} {frac:6.3f} {method:22s} tets {row['tets']:8d} points {row['points']:7d} "
        f"dihedral [{row['min_dihedral']:6.2f}, {row['max_dihedral']:6.2f}] "
        f"<10 {100 * row['under_10']:5.2f}% hausd {row['hausdorff']:.3g} "
        f"dV {100 * row['volume_error']:.3f}% {seconds:7.2f} s {extra}",
        flush=True,
    )


def _mmg3d_file_route(exe, tets, hmin, hmax, hgrad, hausd):
    """mmg3d run as a separate program on a Medit file; its tetrahedra back."""
    with tempfile.TemporaryDirectory() as d:
        src, dst = os.path.join(d, "in.mesh"), os.path.join(d, "out.mesh")
        mp.write(
            src, mp.Mesh(tets.points, [b for b in tets.cells if b.type == "tetra"])
        )
        args = [exe, "-in", src, "-out", dst, "-hmin", repr(hmin), "-hmax", repr(hmax)]
        args += ["-hgrad", repr(hgrad), "-hausd", repr(hausd), "-v", "-1"]
        subprocess.run(args, check=True, capture_output=True, timeout=3600)
        out = mp.read(dst)
        return mp.Mesh(out.points, [b for b in out.cells if b.type == "tetra"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--csv", default=os.path.join(REPO, "benchmark", "adapt_volume.csv")
    )
    ap.add_argument("--hgrad", type=float, default=1.3)
    ap.add_argument(
        "--hausd", type=float, default=0.0, help="0: adapt_volume's default"
    )
    ap.add_argument("--fetch", action="store_true", help="download the classic scans")
    ap.add_argument("--cases", default="sphere,box,cylinder,bunny")
    ap.add_argument(
        "--mmg3d", default=None, help="an mmg3d program, for the file route"
    )
    ap.add_argument("--boundary-sizing", default="input", choices=("input", "geometry"))
    args = ap.parse_args(argv)
    names = args.cases.split(",") + (["buddha", "dragon"] if args.fetch else [])
    have_mmg = getattr(mp, "has_mmg", lambda: False)()
    mmg3d = args.mmg3d or shutil.which("mmg3d_O3") or shutil.which("mmg3d")
    backends = [b for b in mp.tetrahedralize_backends() if b in ("tetgen", "gmsh")]
    rows = []
    for case in names:
        factory, fracs = CASES[case]
        skin = factory()
        if skin is None:
            continue
        diag = _diagonal(skin)
        for frac in fracs:
            min_size = frac * diag
            t0 = time.perf_counter()
            tets, rep = mp.octree_mesh(skin, min_size=min_size, return_report=True)
            t_oct = time.perf_counter() - t0
            q = mp.attach_quality(tets)
            volume_ref = float(
                np.concatenate(
                    [np.asarray(a) for a in q.cell_data["quality:volume"]]
                ).sum()
            )
            _run(rows, case, frac, "octree_mesh", tets, skin, volume_ref, t_oct)
            hmax = rep["max_size"]
            t0 = time.perf_counter()
            out, arep = mp.adapt_volume(
                tets,
                hmin=min_size,
                hmax=hmax,
                hgrad=args.hgrad,
                hausd=args.hausd,
                boundary_sizing=args.boundary_sizing,
                return_report=True,
            )
            t_ad = time.perf_counter() - t0
            _run(
                rows,
                case,
                frac,
                f"adapt_volume {args.boundary_sizing[:4]} {args.hgrad:g}",
                out,
                skin,
                volume_ref,
                t_ad,
                f"hausd {arep['hausd']:.3g} it {arep['iterations']}",
            )
            if have_mmg:
                t0 = time.perf_counter()
                try:
                    mm = mp.adapt(
                        tets,
                        hmin=min_size,
                        hmax=hmax,
                        hgrad=args.hgrad,
                        hausd=arep["hausd"],
                    )
                    _run(
                        rows,
                        case,
                        frac,
                        "mmg3d",
                        mm,
                        skin,
                        volume_ref,
                        time.perf_counter() - t0,
                    )
                except Exception as exc:  # pragma: no cover - optional backend
                    print(f"{case}: mmg3d failed: {exc}")
            elif mmg3d:
                t0 = time.perf_counter()
                try:
                    mm = _mmg3d_file_route(
                        mmg3d, tets, min_size, hmax, args.hgrad, arep["hausd"]
                    )
                    _run(
                        rows,
                        case,
                        frac,
                        "mmg3d (file route)",
                        mm,
                        skin,
                        volume_ref,
                        time.perf_counter() - t0,
                    )
                except Exception as exc:  # pragma: no cover - optional program
                    print(f"{case}: mmg3d failed: {exc}")
            for backend in backends:
                t0 = time.perf_counter()
                try:
                    tg = mp.tetrahedralize(skin, backend=backend, max_size=hmax)
                    _run(
                        rows,
                        case,
                        frac,
                        backend,
                        tg,
                        skin,
                        volume_ref,
                        time.perf_counter() - t0,
                    )
                except Exception as exc:  # pragma: no cover - optional program
                    print(f"{case}: {backend} failed: {exc}")
    if not have_mmg and not mmg3d:
        print("MMG: neither the add-on nor an mmg3d program; its baseline skipped")
    if not backends:
        print("tetgen / gmsh: not installed; their baselines skipped")
    if rows and args.csv:
        os.makedirs(os.path.dirname(args.csv), exist_ok=True)
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
