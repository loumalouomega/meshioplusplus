#  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
# meshio++ — MIT License (see LICENSE). Main authors: Vicente Mataix Ferrandiz
"""Targeted micro-benchmark for the container swaps in the C++ hot paths.

The headline benchmark suite (``bench.py`` / ``01_benchmark.ipynb``) round-trips
formats whose C++ readers/writers only touch the small ``point_data``/
``cell_data`` maps, so it is insensitive to the ``std::map`` -> container changes.
This script instead exercises a genuine per-vertex hot path: the **WKT** reader's
exact-value point dedup, which was ``std::map<std::vector<double>, int64>`` and is
now ``std::unordered_map`` with a custom coordinate hash.

It builds a triangle grid with heavy vertex sharing, writes it to WKT (which emits
every triangle's three vertices explicitly), then times the read — where the
dedup collapses ~3x as many vertex occurrences as there are unique points.

Caveat: WKT read is dominated by text parsing (``strtod`` over tens of MB), so the
dedup container is only a fraction of the wall time and the end-to-end effect of
the container swap is small relative to run-to-run system jitter. To A/B fairly,
build ``main`` and this branch and compare the reported ``min`` (more robust to
jitter than the median) across several invocations of each.

Usage:  python bench_hotpath.py [grid_n]   (default grid_n = 400)
"""

from __future__ import annotations

import os
import statistics
import sys
import tempfile
import time

import numpy as np

import meshioplusplus as pp
from meshioplusplus import Mesh


def _timeit(fn, repeats=15, warmup=2):
    for _ in range(warmup):
        fn()
    ts = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return min(ts), statistics.median(ts)


def triangle_grid(n: int):
    """An ``n x n`` grid of points meshed into ``2*(n-1)**2`` triangles.

    Vertices are shared across triangles, so a WKT round-trip forces the reader
    to dedup ~``6*(n-1)**2`` vertex occurrences down to ``n*n`` unique points.
    """
    lin = np.linspace(0.0, 1.0, n)
    x, y = np.meshgrid(lin, lin, indexing="ij")
    points = np.column_stack([x.ravel(), y.ravel(), np.zeros(n * n)]).astype(np.float64)

    def pid(i, j):
        return i * n + j

    i, j = np.meshgrid(np.arange(n - 1), np.arange(n - 1), indexing="ij")
    i, j = i.ravel(), j.ravel()
    lower = np.stack([pid(i, j), pid(i + 1, j), pid(i + 1, j + 1)], axis=1)
    upper = np.stack([pid(i, j), pid(i + 1, j + 1), pid(i, j + 1)], axis=1)
    tris = np.concatenate([lower, upper], axis=0).astype(np.int64)
    return points, [("triangle", tris)]


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    points, cells = triangle_grid(n)
    ntris = cells[0][1].shape[0]
    print(f"meshio++ {pp.__version__}  |  backend={pp._core.__parallel_backend__}")
    print(
        f"grid n={n}: {len(points):,} unique points, {ntris:,} triangles "
        f"(~{3 * ntris:,} vertex occurrences to dedup on read)"
    )

    mesh = Mesh(points, cells)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "grid.wkt")
        pp.write(path, mesh, file_format="wkt")
        nbytes = os.path.getsize(path)
        read_min, read_med = _timeit(lambda: pp.read(path, file_format="wkt"))
        # sanity: read back the right number of unique points
        got = pp.read(path, file_format="wkt")
        assert got.points.shape[0] == len(points), (got.points.shape[0], len(points))

    print(f"wkt file: {nbytes / 1e6:.1f} MB")
    print(
        f"wkt read (dedup hot path): min={read_min * 1e3:.1f} ms  "
        f"median={read_med * 1e3:.1f} ms  (15 reps)"
    )


if __name__ == "__main__":
    main()
