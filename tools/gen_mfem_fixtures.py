#!/usr/bin/env python3
"""Regenerate the MFEM fixtures under ``tests/python/meshes/mfem/``.

Needs PyMFEM (``pip install mfem``, BSD-3-Clause) and a checkout of MFEM's
``data/`` directory (BSD-3-Clause, https://github.com/mfem/mfem):

    python tools/gen_mfem_fixtures.py /path/to/mfem/data

It does three things, none of which uses meshio++:

1. Copies a few of MFEM's own sample meshes: ``star-q2`` (order-2 quads),
   ``escher-p2`` (order-2 tets), ``fichera-q2`` (order-2 hexes),
   ``fichera-mixed-p2`` (order-2 tets, hexes and prisms), ``compass`` (v1.3
   attribute sets), ``tinyzoo-3d`` (one cell of each 3-D geometry),
   ``periodic-square`` (``L2_T1`` discontinuous nodes), ``escher-p3`` (cubic
   nodes) and ``amr-quad`` (a non-conforming mesh, which must be refused).
2. Projects fields onto some of them with MFEM and saves the grid functions:
   ``star-q2.u.gf`` (H1 P2 scalar), ``star-q2.v.gf`` (H1 P2 2-vector, byNODES),
   ``fichera-q2.w.gf`` (H1 P2 3-vector, byVDIM), ``compass.t.gf`` (H1 P1),
   ``compass.q.gf`` (H1 P2 on a linear mesh) and ``compass.e.gf`` (L2 P0).
3. Freezes MFEM's own evaluation in ``reference.npz``. For every mesh and for
   every vertex, edge, quad face and element of it, keyed by its sorted global
   vertex ids, it stores the point MFEM's element transformation gives at the
   entity's reference centre, and the value of each grid function there. A test
   then needs no MFEM: the node of a meshio++ cell that sits between those
   vertices must be at that point and carry that value.
"""

import pathlib
import shutil
import sys

import numpy as np

try:
    import mfem.ser as mfem
except ImportError:  # pragma: no cover - a maintainer tool
    sys.exit("gen_mfem_fixtures.py needs PyMFEM: pip install mfem")

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "mfem"
)
SAMPLES = [
    "star-q2",
    "escher-p2",
    "fichera-q2",
    "fichera-mixed-p2",
    "compass",
    "tinyzoo-3d",
    "periodic-square",
    "escher-p3",
    "amr-quad",
]
# The meshes the reference covers (the conforming, continuous ones).
REFERENCE = [
    "star-q2",
    "escher-p2",
    "fichera-q2",
    "fichera-mixed-p2",
    "compass",
    "tinyzoo-3d",
]

# MFEM's local edges and faces per geometry (fem/geom.cpp), MFEM vertex order.
EDGES = {
    2: [(0, 1), (1, 2), (2, 0)],
    3: [(0, 1), (1, 2), (2, 3), (3, 0)],
    4: [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
    5: [(0, 1), (1, 2), (3, 2), (0, 3), (4, 5), (5, 6), (7, 6), (4, 7), (0, 4), (1, 5)]
    + [(2, 6), (3, 7)],
    6: [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    7: [(0, 1), (1, 2), (3, 2), (0, 3), (0, 4), (1, 4), (2, 4), (3, 4)],
}
QUAD_FACES = {
    5: [
        (3, 2, 1, 0),
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
        (4, 5, 6, 7),
    ],
    6: [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    7: [(3, 2, 1, 0)],
}


class Scalar(mfem.PyCoefficient):
    def __init__(self, f):
        super().__init__()
        self.f = f

    def EvalValue(self, x):
        return self.f(x)


class Vector(mfem.VectorPyCoefficient):
    def __init__(self, dim, f):
        super().__init__(dim)
        self.f = f

    def EvalValue(self, x):
        return self.f(x)


def project(mesh, order, vdim, ordering, f, l2=False):
    dim = mesh.Dimension()
    fec = mfem.L2_FECollection(order, dim) if l2 else mfem.H1_FECollection(order, dim)
    fes = mfem.FiniteElementSpace(mesh, fec, vdim, ordering)
    gf = mfem.GridFunction(fes)
    gf.ProjectCoefficient(Scalar(f) if vdim == 1 else Vector(vdim, f))
    # keep the objects alive with the grid function
    gf._keep = (fec, fes)
    return gf


def fields(name, mesh):
    by_nodes, by_vdim = mfem.Ordering.byNODES, mfem.Ordering.byVDIM
    if name == "star-q2":
        return {
            "u": project(mesh, 2, 1, by_nodes, lambda x: x[0] * x[0] + 2 * x[1]),
            "v": project(mesh, 2, 2, by_nodes, lambda x: (x[0] - x[1], x[0] * x[1])),
        }
    if name == "fichera-q2":
        return {
            "w": project(
                mesh, 2, 3, by_vdim, lambda x: (x[0] * x[1], x[1] + x[2], x[2] * x[2])
            )
        }
    if name == "compass":
        return {
            "t": project(mesh, 1, 1, by_nodes, lambda x: 3 * x[0] - x[1] + 1),
            "q": project(mesh, 2, 1, by_nodes, lambda x: x[0] * x[1] + x[0] * x[0]),
            "e": project(mesh, 0, 1, by_nodes, lambda x: x[0] + 10 * x[1], l2=True),
        }
    return {}


def entities(mesh, i):
    """(sorted global vertex key, reference point) of every vertex, edge, quad face
    and the element itself."""
    geom = mesh.GetElementBaseGeometry(i)
    v = list(mesh.GetElementVertices(i))
    rule = mfem.Geometries.GetVertices(geom)
    ref = np.array(
        [
            [rule.IntPoint(k).x, rule.IntPoint(k).y, rule.IntPoint(k).z]
            for k in range(len(v))
        ]
    )
    subsets = (
        [(k,) for k in range(len(v))] + EDGES.get(geom, []) + QUAD_FACES.get(geom, [])
    )
    subsets.append(tuple(range(len(v))))
    for s in subsets:
        yield tuple(sorted(v[k] for k in s)), ref[list(s)].mean(0)


def reference(name, mesh, gfs):
    keys, xyz, values = {}, [], {g: [] for g in gfs}
    sdim = mesh.SpaceDimension()
    for i in range(mesh.GetNE()):
        T = mesh.GetElementTransformation(i)
        for key, point in entities(mesh, i):
            if key in keys:
                continue
            ip = mfem.IntegrationPoint()
            ip.Set3(*point)
            T.SetIntPoint(ip)
            keys[key] = len(xyz)
            x = np.zeros(3)
            x[:sdim] = np.array(T.Transform(ip))
            xyz.append(x)
            for g, gf in gfs.items():
                vdim = gf.FESpace().GetVDim()
                if vdim == 1:
                    values[g].append([gf.GetValue(T, ip)])
                else:
                    out = mfem.Vector()
                    gf.GetVectorValue(T, ip, out)
                    values[g].append(list(out.GetDataArray()))
    padded = np.full((len(keys), 8), -1, dtype=np.int64)
    for key, k in keys.items():
        padded[k, : len(key)] = key
    arrays = {f"{name}/keys": padded, f"{name}/xyz": np.array(xyz)}
    for g, vals in values.items():
        arrays[f"{name}/{g}"] = np.array(vals)
    return arrays


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    data = pathlib.Path(sys.argv[1])
    OUT.mkdir(parents=True, exist_ok=True)
    for name in SAMPLES:
        shutil.copyfile(data / f"{name}.mesh", OUT / f"{name}.mesh")
    arrays = {}
    for name in REFERENCE:
        mesh = mfem.Mesh(str(OUT / f"{name}.mesh"), 1, 1)
        gfs = fields(name, mesh)
        for g, gf in gfs.items():
            gf.Save(str(OUT / f"{name}.{g}.gf"))
        arrays.update(reference(name, mesh, gfs))
    np.savez_compressed(OUT / "reference.npz", **arrays)


if __name__ == "__main__":
    main()
