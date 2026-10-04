#  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
# meshio++ — MIT License (see LICENSE). Main authors: Vicente Mataix Ferrandiz
"""Per-call cost of the Python <-> C++ boundary (roadmap 3.4.3).

Three cases, each a ``Mesh`` passed into the core and the result built back, so
both directions of ``np_conversions.hpp`` are timed and no file I/O is involved:

* ``tiny``: a two-tetrahedron mesh through ``clean`` 20 000 times -- the fixed
  cost of ``py_to_mesh`` + ``mesh_to_py`` (class lookup, byte-order checks);
* ``polygons``: a ``polygon`` block of N triangles-as-ragged rows (N = 200 000);
* ``polyhedra``: a ``polyhedron`` block of N hexahedra (N = 50 000).

Usage::

    python benchmark/bench_boundary.py [--repeats 5]

Run it on the base and on the changed build and compare the ``min`` column,
which is more robust to jitter than the median.
"""

from __future__ import annotations

import argparse
import statistics
import time

import numpy as np

import meshioplusplus as pp
from meshioplusplus import Mesh


def _timeit(fn, repeats, warmup=1):
    for _ in range(warmup):
        fn()
    ts = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return min(ts), statistics.median(ts)


def tiny_mesh():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], float)
    return Mesh(pts, [("tetra", np.array([[0, 1, 2, 3], [1, 2, 3, 4]], np.int64))])


def polygon_mesh(n):
    # n polygons: quads and triangles alternating, on a strip of shared points
    m = n + 2
    pts = np.zeros((2 * m, 3))
    pts[:m, 0] = np.arange(m)
    pts[m:, 0] = np.arange(m)
    pts[m:, 1] = 1.0
    rows = []
    for i in range(n):
        if i % 2:
            rows.append(np.array([i, i + 1, m + i + 1, m + i], np.int64))
        else:
            rows.append(np.array([i, i + 1, m + i], np.int64))
    return Mesh(pts, [("polygon", rows)])


def polyhedron_mesh(n):
    # n disjoint unit hexahedra, six quad faces each
    corner = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        float,
    )
    pts = (
        corner[None] + 2.0 * np.arange(n)[:, None, None] * np.array([1.0, 0, 0])
    ).reshape(-1, 3)
    faces = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ]
    cells = [[np.array(f, np.int64) + 8 * i for f in faces] for i in range(n)]
    return Mesh(pts, [("polyhedron", cells)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--polygons", type=int, default=200_000)
    ap.add_argument("--polyhedra", type=int, default=50_000)
    args = ap.parse_args()

    tiny = tiny_mesh()
    cases = [
        ("tiny x20000", lambda: [pp.clean(tiny) for _ in range(20_000)]),
    ]
    poly = polygon_mesh(args.polygons)
    cases.append((f"polygons n={args.polygons}", lambda: pp.clean(poly)))
    polyh = polyhedron_mesh(args.polyhedra)
    cases.append((f"polyhedra n={args.polyhedra}", lambda: pp.clean(polyh)))

    print(f"{'case':<26}{'min [s]':>10}{'median [s]':>12}")
    for name, fn in cases:
        lo, med = _timeit(fn, args.repeats)
        print(f"{name:<26}{lo:>10.4f}{med:>12.4f}")


if __name__ == "__main__":
    main()
