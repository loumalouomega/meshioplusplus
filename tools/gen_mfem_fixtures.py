#!/usr/bin/env python3
"""Regenerate the MFEM fixtures under ``tests/python/meshes/mfem/``.

Needs PyMFEM (``pip install mfem``, BSD-3-Clause) and a checkout of MFEM's
``data/`` directory (BSD-3-Clause, https://github.com/mfem/mfem):

    python tools/gen_mfem_fixtures.py /path/to/mfem/data

It does four things, none of which uses meshio++:

1. Copies a few of MFEM's own sample meshes: ``star-q2`` (order-2 quads),
   ``escher-p2`` (order-2 tets), ``fichera-q2`` (order-2 hexes),
   ``fichera-mixed-p2`` (order-2 tets, hexes and prisms), ``compass`` (v1.3
   attribute sets), ``tinyzoo-3d`` (one cell of each 3-D geometry),
   ``periodic-square`` (``L2_T1`` discontinuous nodes), ``escher-p3`` and
   ``fichera-q3`` (legacy ``Cubic`` tets and hexes), ``toroid-wedge`` (``H1``
   order-3 prisms), ``rt-2d-p4-tri`` (order-4 triangles) and ``amr-quad`` (an
   ``MFEM NC mesh``).
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
4. For the arbitrary-order meshes -- the MFEM samples above of order 3 and 4,
   and ``curved-*`` meshes it makes by curving MFEM's one-element reference
   meshes to order 3-5 (Gauss-Lobatto, and equispaced ``H1@U`` for the quad)
   and warping them with a smooth non-polynomial map -- it saves an order-3 to
   5 ``H1`` field ``<name>.u.gf`` and freezes MFEM's own high-order output
   (``ParaViewDataCollection`` with ``SetHighOrderOutput``, VTK Lagrange cells)
   in ``reference_lagrange.npz``: per case and cell type, each cell's nodes in
   VTK order and ``u`` at them. For ``amr-quad`` it stores the attribute and
   sorted corner coordinates of every leaf and boundary element MFEM builds.
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
    "fichera-q3",
    "toroid-wedge",
    "rt-2d-p4-tri",
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


# --- arbitrary order --------------------------------------------------------------

# (name, source mesh, curvature order or None, basis, uniform refinements,
# order of the u field or None): MFEM's own order-3/4 samples as they are, and
# MFEM's one-element reference meshes curved and warped. Loading a mesh of
# triangles or tetrahedra re-marks them (reorders their vertices), so a field
# saved against it only matches a mesh printed after that load: the escher-p3
# and rt-2d-p4-tri samples stay MFEM's files and get no field.
HIGH_ORDER = [
    ("escher-p3", "escher-p3", None, None, 0, None),
    ("fichera-q3", "fichera-q3", None, None, 0, 3),
    ("toroid-wedge", "toroid-wedge", None, None, 0, 3),
    ("rt-2d-p4-tri", "rt-2d-p4-tri", None, None, 0, None),
    ("curved-tet-p4", "ref-tetrahedron", 4, "G", 1, 2),
    ("curved-hex-p3", "ref-cube", 3, "G", 1, 5),
    ("curved-prism-p4", "ref-prism", 4, "G", 1, 3),
    ("curved-tri-p5", "ref-triangle", 5, "G", 1, 3),
    ("curved-quad-p3u", "ref-square", 3, "U", 1, 3),
    ("curved-segment-p3", "ref-segment", 3, "G", 2, 3),
]


def _warp(x):
    """A smooth, non-polynomial map, so every order is exercised."""
    y = np.array(x, dtype=float)
    d = len(y)
    y[0] = x[0] + 0.1 * np.sin(1.3 * (x[1] if d > 1 else x[0]))
    if d > 1:
        y[1] = x[1] + 0.07 * np.cos(0.9 * x[0])
    if d > 2:
        y[2] = x[2] + 0.05 * np.sin(x[0] + x[1])
    return y


def _field(x):
    return float(np.sin(x[0]) + x[-1] ** 3 + 0.5 * np.cos(sum(x)))


def _vtu_cells(path):
    """(type code, node coordinates, u) of every cell of an ASCII VTU."""
    import xml.etree.ElementTree as ET

    root = ET.parse(path).getroot()
    arrays = {}
    for da in root.iter("DataArray"):
        values = np.array(da.text.split(), dtype=float)
        arrays[da.get("Name") or "Points"] = (
            values,
            int(da.get("NumberOfComponents") or 1),
        )
    pts = arrays["Points"][0].reshape(-1, 3)
    conn = arrays["connectivity"][0].astype(int)
    offsets = arrays["offsets"][0].astype(int)
    types = arrays["types"][0].astype(int)
    u = arrays["u"][0] if "u" in arrays else None
    out, start = [], 0
    for t, end in zip(types, offsets):
        ids = conn[start:end]
        out.append((int(t), pts[ids], None if u is None else u[ids]))
        start = end
    return out


def high_order(data, name, source, curve, basis, refine, u_order, arrays):
    import os
    import tempfile

    mesh = mfem.Mesh(str(data / f"{source}.mesh"), 1, 1)
    if curve:
        for _ in range(refine):
            mesh.UniformRefinement()
        btype = (
            mfem.BasisType.ClosedUniform
            if basis == "U"
            else mfem.BasisType.GaussLobatto
        )
        fec = mfem.H1_FECollection(curve, mesh.Dimension(), btype)
        fes = mfem.FiniteElementSpace(
            mesh, fec, mesh.SpaceDimension(), mfem.Ordering.byVDIM
        )
        mesh.SetNodalFESpace(fes)

        class Warp(mfem.VectorPyCoefficient):
            def EvalValue(self, x):
                return _warp(x)

        mesh.GetNodes().ProjectCoefficient(Warp(mesh.SpaceDimension()))
    path = OUT / f"{name}.mesh"
    if curve:
        mesh.Print(str(path), 17)
        # Loading re-marks tetrahedra (reorders their vertices): save the mesh
        # as it then stands, so the file and the field agree.
        mesh = mfem.Mesh(str(path), 1, 1)
        mesh.Print(str(path), 17)
    else:
        mesh = mfem.Mesh(str(path), 1, 1)

    class Field(mfem.PyCoefficient):
        def EvalValue(self, x):
            return _field(x)

    order = mesh.GetNodes().FESpace().GetMaxElementOrder()
    fec_u = mfem.H1_FECollection(u_order or 1, mesh.Dimension())
    fes_u = mfem.FiniteElementSpace(mesh, fec_u)
    u = mfem.GridFunction(fes_u)
    u.ProjectCoefficient(Field())
    if u_order:
        u.Save(str(OUT / f"{name}.u.gf"), 17)
        order = max(order, u_order)
    with tempfile.TemporaryDirectory() as tmp:
        dc = mfem.ParaViewDataCollection(name, mesh)
        dc.SetPrefixPath(tmp)
        dc.SetLevelsOfDetail(order)
        dc.SetHighOrderOutput(True)
        dc.SetDataFormat(mfem.VTKFormat_ASCII)
        dc.SetPrecision(17)
        if u_order:
            dc.RegisterField("u", u)
        dc.Save()
        vtu = next(
            os.path.join(r, f)
            for r, _, fs in os.walk(tmp)
            for f in fs
            if f.endswith(".vtu")
        )
        cells = _vtu_cells(vtu)
    for code in sorted({c[0] for c in cells}):
        group = [c for c in cells if c[0] == code]
        dim = mesh.SpaceDimension()
        arrays[f"{name}:{code}:points"] = np.array([c[1][:, :dim] for c in group])
        if u_order:
            arrays[f"{name}:{code}:u"] = np.array([c[2] for c in group])
    print(name, mesh.GetNE(), "cells, VTK order", order)


def non_conforming(name, arrays):
    mesh = mfem.Mesh(str(OUT / f"{name}.mesh"), 1, 1)
    sdim = mesh.SpaceDimension()

    def cells(n, vertices, attribute):
        corners, attrs = [], []
        for i in range(n):
            pts = sorted(tuple(mesh.GetVertexArray(v)[:sdim]) for v in vertices(i))
            corners.append(pts)
            attrs.append(attribute(i))
        return np.array(corners), np.array(attrs)

    for part, n, verts, attr in (
        ("elements", mesh.GetNE(), mesh.GetElementVertices, mesh.GetAttribute),
        ("boundary", mesh.GetNBE(), mesh.GetBdrElementVertices, mesh.GetBdrAttribute),
    ):
        corners, attrs = cells(n, verts, attr)
        arrays[f"{name}:{part}:corners"] = corners
        arrays[f"{name}:{part}:attributes"] = attrs


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
    lagrange = {}
    for case in HIGH_ORDER:
        high_order(data, *case, lagrange)
    non_conforming("amr-quad", lagrange)
    np.savez_compressed(OUT / "reference_lagrange.npz", **lagrange)


if __name__ == "__main__":
    main()
