#  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
# meshio++ — MIT License (see LICENSE). Main authors: Vicente Mataix Ferrandiz
"""What the browser viewer's zlib-compressed VTP costs to produce and consume (roadmap 3.4.5.1).

``renderPipeline`` in the viewer's worker writes every operation result as a
zlib-compressed binary VTP, and vtk.js base64-decodes and inflates it on the
main thread. This script measures, on a triangulated surface carrying a point
scalar, a cell vector and an int64 provenance array, for each of the three
encodings the roadmap discusses (zlib, uncompressed with the default 4-byte
size header, uncompressed with the 8-byte header vtk.js could read):

* ``write``: the cost of ``meshioplusplus.vtp.write`` (producer side);
* ``b64``: base64-decoding every ``<DataArray>`` of the file;
* ``inflate``: ``zlib`` inflating the decoded blocks (zero for the uncompressed
  encodings), the part the roadmap proposes to skip;
* the file size, which is what crosses ``postMessage``.

Usage::

    python benchmark/bench_vtp_inflate.py [--sizes 100,300,700] [--repeats 5] [--csv out.csv]

``--sizes`` are grid sides ``n`` of an ``n x n`` quad grid split into
``2 n^2`` triangles. Compare the ``min`` columns; they are more robust to
jitter than the median. Use ``benchmark/vtp_inflate.js`` for vtk.js's own
parse time on the same files (``--keep DIR`` writes them).
"""

from __future__ import annotations

import argparse
import base64
import csv
import re
import statistics
import struct
import sys
import tempfile
import time
import zlib
from pathlib import Path

import numpy as np

import meshioplusplus as pp
from meshioplusplus import Mesh

_ARRAY = re.compile(r"<DataArray\b[^>]*>(.*?)</DataArray>", re.S)


def _timeit(fn, repeats, warmup=1):
    for _ in range(warmup):
        fn()
    ts = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return min(ts), statistics.median(ts)


def surface_mesh(n):
    """An ``n x n`` quad grid split into triangles, with the arrays the viewer keeps."""
    xs, ys = np.meshgrid(np.linspace(0, 1, n + 1), np.linspace(0, 1, n + 1))
    pts = np.column_stack([xs.ravel(), ys.ravel(), np.sin(3 * xs.ravel()) * ys.ravel()])
    idx = np.arange((n + 1) * (n + 1)).reshape(n + 1, n + 1)
    a, b = idx[:-1, :-1].ravel(), idx[:-1, 1:].ravel()
    c, d = idx[1:, 1:].ravel(), idx[1:, :-1].ravel()
    tris = np.concatenate([np.column_stack([a, b, c]), np.column_stack([a, c, d])]).astype(np.int64)
    nc = len(tris)
    rng = np.random.default_rng(0)
    return Mesh(
        pts,
        [("triangle", tris)],
        point_data={"scalar": pts[:, 2].copy()},
        cell_data={
            "vector": [rng.standard_normal((nc, 3))],
            "surface:parent_cell": [np.arange(nc, dtype=np.int64)],
        },
    )


def _payloads(path):
    """The base64 text of each ``<DataArray>`` of a VTP file."""
    text = Path(path).read_text()
    return [m.group(1).strip() for m in _ARRAY.finditer(text) if m.group(1).strip()]


def _decode(payloads, header_size):
    """Base64-decode every array and split its header; returns the raw pieces."""
    out = []
    for p in payloads:
        raw = base64.b64decode(p)
        out.append(raw)
    return out


def _inflate(raws, header_size):
    """Inflate the blocks of every zlib array (vtk.js's second step)."""
    total = 0
    for raw in raws:
        hdr_t = "<I" if header_size == 4 else "<Q"
        w = struct.calcsize(hdr_t)
        nblocks, _bs, _last = struct.unpack_from("<" + hdr_t[1] * 3, raw, 0)
        sizes = struct.unpack_from("<" + hdr_t[1] * nblocks, raw, 3 * w)
        pos = (3 + nblocks) * w
        for s in sizes:
            total += len(zlib.decompress(raw[pos : pos + s]))
            pos += s
    return total


def _variants():
    return [
        # label, writer kwargs, header bytes, compressed?
        ("zlib", dict(compression="zlib"), 4, True),
        ("raw4", dict(compression=None), 4, False),
        ("raw8", dict(compression=None, header_type="UInt64"), 8, False),
    ]


def measure(n, repeats, outdir):
    mesh = surface_mesh(n)
    rows = []
    for label, kwargs, hsz, compressed in _variants():
        path = Path(outdir) / f"surf_{n}_{label}.vtp"
        write_min, write_med = _timeit(lambda: pp.vtp.write(path, mesh, **kwargs), repeats)
        payloads = _payloads(path)
        b64_min, b64_med = _timeit(lambda: _decode(payloads, hsz), repeats)
        raws = _decode(payloads, hsz)
        if compressed:
            inf_min, inf_med = _timeit(lambda: _inflate(raws, hsz), repeats)
        else:
            inf_min = inf_med = 0.0
        rows.append(
            {
                "n": n,
                "cells": 2 * n * n,
                "encoding": label,
                "bytes": path.stat().st_size,
                "write_min_ms": write_min * 1e3,
                "b64_min_ms": b64_min * 1e3,
                "inflate_min_ms": inf_min * 1e3,
                "write_med_ms": write_med * 1e3,
                "b64_med_ms": b64_med * 1e3,
                "inflate_med_ms": inf_med * 1e3,
            }
        )
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sizes", default="100,300,700", help="comma-separated grid sides")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--csv", help="also write the rows to this CSV file")
    ap.add_argument("--keep", help="write the VTP files here and keep them (for vtp_inflate.js)")
    args = ap.parse_args(argv)

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        outdir = Path(args.keep) if args.keep else Path(tmp)
        outdir.mkdir(parents=True, exist_ok=True)
        for n in [int(s) for s in args.sizes.split(",")]:
            rows.extend(measure(n, args.repeats, outdir))

    cols = list(rows[0])
    print("  ".join(f"{c:>14}" for c in cols))
    for r in rows:
        print("  ".join(f"{r[c]:>14.3f}" if isinstance(r[c], float) else f"{r[c]:>14}" for c in cols))
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
