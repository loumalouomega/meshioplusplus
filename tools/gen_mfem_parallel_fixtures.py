#!/usr/bin/env python3
"""Regenerate the parallel MFEM fixtures under ``tests/python/meshes/mfem/parallel/``.

Builds ``tools/mfem_parallel_driver.cpp`` against MFEM compiled with MPI (and
hypre), runs it with ``mpirun`` and freezes MFEM's own high-order output of each
case in ``reference.npz``, so the tests need no MFEM:

    MFEM_BUILD=/path/to/mfem/build-par MFEM_SRC=/path/to/mfem \\
    HYPRE_DIR=/path/to/hypre-install \\
        python tools/gen_mfem_parallel_fixtures.py /path/to/mfem/data

MFEM (BSD-3-Clause) was configured with ``-DMFEM_USE_MPI=YES
-DMFEM_USE_METIS=NO`` (the driver partitions itself). Cases:

* ``star-p2``: MFEM's ``star.mesh`` refined once, curved to order 2 and warped,
  over 4 ranks, with an order-2 field ``u``.
* ``beam-tet``: MFEM's ``beam-tet.mesh`` refined once, linear, over 3 ranks, with
  an order-1 field ``u``.
* ``star-nc``: MFEM's ``star.mesh`` refined once, made non-conforming and refined
  twice more in bands, over 3 ranks, with an order-2 field ``u``: its
  ``.pmesh.NNNNNN`` files are ``MFEM NC mesh`` files (ParPrint of a
  non-conforming mesh), each rank's refinement tree with its ghosts.

Each has both layouts: ``<case>.mesh.NNNNNN`` (``ParMesh::Save``) and
``<case>.pmesh.NNNNNN`` (``ParMesh::ParPrint``, with communication groups). For
each case ``reference.npz`` holds, per VTK cell type, every rank's cells' nodes
(``<case>:<type>:points``) and ``u`` at them (``<case>:<type>:u``).
"""

import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent / "tests" / "python" / "meshes" / "mfem" / "parallel"
CASES = [
    ("star-p2", "star", 2, 2, 4, False),
    ("beam-tet", "beam-tet", 0, 1, 3, False),
    ("star-nc", "star", 0, 2, 3, True),
]


def _vtu_cells(path):
    root = ET.parse(path).getroot()
    arrays = {}
    for da in root.iter("DataArray"):
        arrays[da.get("Name") or "Points"] = np.array(da.text.split(), dtype=float)
    pts = arrays["Points"].reshape(-1, 3)
    conn = arrays["connectivity"].astype(int)
    offsets = arrays["offsets"].astype(int)
    out, start = [], 0
    for t, end in zip(arrays["types"].astype(int), offsets):
        ids = conn[start:end]
        out.append((int(t), pts[ids], arrays["u"][ids]))
        start = end
    return out


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    data = pathlib.Path(sys.argv[1])
    build, src, hypre = (os.environ[k] for k in ("MFEM_BUILD", "MFEM_SRC", "HYPRE_DIR"))
    OUT.mkdir(parents=True, exist_ok=True)
    arrays = {}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        driver = tmp / "driver"
        subprocess.run(
            [
                "mpicxx",
                "-O1",
                "-std=c++17",
                f"-I{build}",
                f"-I{src}",
                f"-I{hypre}/include",
                str(HERE / "mfem_parallel_driver.cpp"),
                f"{build}/libmfem.a",
                f"{hypre}/lib/libHYPRE.a",
                "-lm",
                "-o",
                str(driver),
            ],
            check=True,
        )
        for name, source, curve, u_order, ranks, nc in CASES:
            prefix = tmp / name
            subprocess.run(
                [
                    "mpirun",
                    "--oversubscribe",
                    "-np",
                    str(ranks),
                    str(driver),
                    str(data / f"{source}.mesh"),
                    str(prefix),
                    str(curve),
                    str(u_order),
                ]
                + (["nc"] if nc else []),
                check=True,
            )
            for r in range(ranks):
                for kind in ("mesh", "pmesh", "u"):
                    shutil.copyfile(
                        f"{prefix}.{kind}.{r:06d}", OUT / f"{name}.{kind}.{r:06d}"
                    )
            cells = []
            for vtu in sorted((tmp / f"{name}_ref").rglob("proc*.vtu")):
                cells += _vtu_cells(vtu)
            for code in sorted({c[0] for c in cells}):
                group = [c for c in cells if c[0] == code]
                arrays[f"{name}:{code}:points"] = np.array([c[1] for c in group])
                arrays[f"{name}:{code}:u"] = np.array([c[2] for c in group])
            print(name, len(cells), "cells over", ranks, "ranks")
    np.savez_compressed(OUT / "reference.npz", **arrays)


if __name__ == "__main__":
    main()
