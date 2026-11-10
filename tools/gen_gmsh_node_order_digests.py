#!/usr/bin/env python3
"""Record Gmsh's own node order of its higher-order elements, as digests.

meshio++ holds the fixed-order types beyond VTK's (``triangle10``,
``quad16``, ``tetra20``, ``wedge40``, ``hexahedron64``, ...) in Gmsh's node
order and writes them as VTK Lagrange cells through
``_lagrange.gmsh_to_vtk_order``, which generates that order. This script asks
Gmsh itself (``pip install gmsh``; it runs Gmsh in this process, nothing of it
is copied) for every element it has and prints one digest per ``shape:order``
of the node list as lattice coordinates, in ``_lagrange.vtk_lattice``'s frame;
``tests/python/test_lagrange_gmsh_order.py`` holds the generator to them.

    python tools/gen_gmsh_node_order_digests.py        # print the table

``getElementProperties`` describes every element but the prism of order 3 and
up ("FaceClosureFull not implemented"); those are read from a meshed prism
raised to the order instead.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np

FAMILY = {
    "line": "Line",
    "triangle": "Triangle",
    "quad": "Quadrangle",
    "tetra": "Tetrahedron",
    "hexahedron": "Hexahedron",
    "wedge": "Prism",
}
#: Which reference coordinates run over [-1, 1] (the rest over [0, 1]).
SYMMETRIC = {
    "line": (True,),
    "triangle": (False, False),
    "quad": (True, True),
    "tetra": (False, False, False),
    "hexahedron": (True, True, True),
    "wedge": (False, False, True),
}
#: The orders the meshio++ types cover, from 3 (below it VTK has its own ids).
ORDERS = {"hexahedron": range(3, 10), "wedge": range(3, 10)}


def digest(nodes):
    """The digest of a node list ``[(i, j, k), ...]``."""
    text = json.dumps([[int(x) for x in n] for n in nodes])
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def _from_properties(gmsh, shape, p):
    kind = gmsh.model.mesh.getElementType(FAMILY[shape], p)
    _, dim, _, n, coords, _ = gmsh.model.mesh.getElementProperties(kind)
    out = []
    for row in np.asarray(coords).reshape(n, dim):
        ijk = [
            round((x + 1) / 2 * p) if sym else round(x * p)
            for x, sym in zip(row, SYMMETRIC[shape])
        ]
        out.append(tuple(int(v) for v in ijk) + (0,) * (3 - dim))
    return out


def _from_meshed_prism(gmsh, p):
    gmsh.model.add(f"prism{p}")
    geo = gmsh.model.geo
    a, b, c = (geo.addPoint(*x) for x in ((0, 0, 0), (1, 0, 0), (0, 1, 0)))
    lines = [geo.addLine(a, b), geo.addLine(b, c), geo.addLine(c, a)]
    for line in lines:
        geo.mesh.setTransfiniteCurve(line, 2)
    surface = geo.addPlaneSurface([geo.addCurveLoop(lines)])
    geo.extrude([(2, surface)], 0, 0, 1, numElements=[1], recombine=True)
    geo.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(p)
    _, _, nodes = gmsh.model.mesh.getElements(3)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    position = dict(zip(tags, np.asarray(coords).reshape(-1, 3)))
    x = np.array([position[t] for t in nodes[0]])
    frame = np.array([x[1] - x[0], x[2] - x[0], x[3] - x[0]]).T
    uvw = np.linalg.solve(frame, (x - x[0]).T).T
    gmsh.model.remove()
    return [tuple(int(round(v * p)) for v in r) for r in uvw]


def gmsh_digests():
    """``{"shape:order": digest}`` from Gmsh."""
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        out = {}
        for shape in FAMILY:
            for p in ORDERS.get(shape, range(3, 11)):
                if shape == "wedge":
                    nodes = _from_meshed_prism(gmsh, p)
                else:
                    nodes = _from_properties(gmsh, shape, p)
                out[f"{shape}:{p}"] = digest(nodes)
        return out
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    for key, value in gmsh_digests().items():
        print(f'    "{key}": "{value}",')
