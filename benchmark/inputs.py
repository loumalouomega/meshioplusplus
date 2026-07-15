"""Benchmark input meshes.

Two inputs:

* ``example_geometry()`` -- the bundled ``example/example.msh`` bracket reduced
  to its portable geometry (surface triangles + solid tetrahedra).
* ``synthetic_tet_grid(n)`` -- a numpy-generated structured tetrahedral mesh
  (a unit cube split into 6 tets per voxel), sized to make read/write timings
  meaningful. Regenerated on demand; nothing is cached to disk.

Both return ``(points, cells)`` where ``cells`` is a list of
``(type, connectivity)`` tuples accepted by *both* ``meshio`` and
``meshioplusplus``.
"""

from __future__ import annotations

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.normpath(os.path.join(HERE, "..", "example", "example.msh"))


def example_geometry():
    """Bracket geometry (triangles + tetrahedra) from example.msh."""
    import meshioplusplus as mp

    src = mp.read(EXAMPLE)
    tri = np.concatenate([cb.data for cb in src.cells if cb.type == "triangle"])
    tet = np.concatenate([cb.data for cb in src.cells if cb.type == "tetra"])
    return src.points, [("triangle", tri), ("tetra", tet)]


# The 6-tetrahedron (Freudenthal) split of a cube, as corner-index triples into
# the 8 voxel corners c0..c7 (c0=(0,0,0), c1=(1,0,0), c2=(0,1,0), c3=(1,1,0),
# c4=(0,0,1), ...): a standard consistent decomposition.
_TET6 = np.array(
    [
        [0, 1, 3, 7],
        [0, 3, 2, 7],
        [0, 2, 6, 7],
        [0, 6, 4, 7],
        [0, 4, 5, 7],
        [0, 5, 1, 7],
    ],
    dtype=np.int64,
)


def synthetic_tet_grid(n: int = 56):
    """A structured tetrahedral mesh of a unit cube, ~``6*(n-1)**3`` tets."""
    lin = np.linspace(0.0, 1.0, n)
    x, y, z = np.meshgrid(lin, lin, lin, indexing="ij")
    points = np.column_stack([x.ravel(), y.ravel(), z.ravel()]).astype(np.float64)

    def pid(i, j, k):
        return (i * n + j) * n + k

    m = n - 1
    i, j, k = np.meshgrid(np.arange(m), np.arange(m), np.arange(m), indexing="ij")
    i, j, k = i.ravel(), j.ravel(), k.ravel()

    # 8 corners of every voxel, ordered c0..c7 to match _TET6.
    corners = np.stack(
        [
            pid(i, j, k),
            pid(i + 1, j, k),
            pid(i, j + 1, k),
            pid(i + 1, j + 1, k),
            pid(i, j, k + 1),
            pid(i + 1, j, k + 1),
            pid(i, j + 1, k + 1),
            pid(i + 1, j + 1, k + 1),
        ],
        axis=1,
    )  # (nvoxel, 8)

    # (nvoxel, 6, 4) -> (nvoxel*6, 4)
    tets = corners[:, _TET6].reshape(-1, 4).astype(np.int64)
    return points, [("tetra", tets)]
