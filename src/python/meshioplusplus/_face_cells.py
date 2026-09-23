"""Volume cells from their faces: the kernel shared by the readers of
face-based formats (OpenFOAM ``polyMesh``, ANSYS Fluent ``.msh``).

A cell is given as its faces, each a list of point ids wound so that its
right-hand normal points out of the cell. Tetrahedra, pyramids, wedges and
hexahedra come back in meshio++'s (VTK) node order with a positive volume
(the orientation is checked geometrically, so a flipped input face does not
invert the cell); anything else stays a polyhedron of its outward faces.
Twin of ``src/cpp/src/formats/face_cells_common.hpp``.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from ._mesh import CellBlock

__all__ = [
    "triple",
    "node_adjacency",
    "match_top",
    "build_tetra",
    "build_pyramid",
    "build_wedge",
    "build_hexahedron",
    "build_polyhedra",
    "reconstruct_cell",
    "polygon_from_edges",
]


def triple(a, b, c) -> float:
    """Scalar triple product a · (b × c)."""
    return float(np.dot(a, np.cross(b, c)))


def node_adjacency(faces) -> dict:
    """Build a node-to-node adjacency dict from a list of faces."""
    adj: dict[int, set] = {}
    for f in faces:
        m = len(f)
        for i in range(m):
            a, b = f[i], f[(i + 1) % m]
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
    return adj


def match_top(bottom, oriented):
    """
    For each base node, find its unique vertical neighbour.
    Returns the ordered top ring, or None if the topology is ambiguous.
    """
    adj = node_adjacency(oriented)
    base = set(bottom)
    top = []
    for b in bottom:
        cand = [x for x in adj[b] if x not in base]
        if len(cand) != 1:
            return None
        top.append(cand[0])
    return top


def build_tetra(oriented, P):
    """Build a tetrahedron connectivity with positive volume orientation."""
    base = oriented[0]
    apex = (set().union(*oriented) - set(base)).pop()
    n = [base[0], base[1], base[2], apex]
    p = [P[i] for i in n]
    if triple(p[1] - p[0], p[2] - p[0], p[3] - p[0]) < 0:
        n = [base[0], base[2], base[1], apex]
    return n


def build_pyramid(oriented, P):
    """Build a pyramid connectivity with positive volume orientation."""
    quad = next(f for f in oriented if len(f) == 4)
    apex = (set().union(*oriented) - set(quad)).pop()
    n = list(quad) + [apex]
    p = [P[i] for i in n]
    if triple(p[1] - p[0], p[3] - p[0], p[4] - p[0]) < 0:
        n = [quad[0], quad[3], quad[2], quad[1], apex]
    return n


def build_wedge(oriented, P):
    """Build a wedge connectivity with positive volume orientation."""
    bottom = next(f for f in oriented if len(f) == 3)
    top = match_top(bottom, oriented)
    if top is None:
        return None
    n = list(bottom) + top
    p = [P[i] for i in n]
    if triple(p[1] - p[0], p[2] - p[0], p[3] - p[0]) < 0:
        n = [bottom[0], bottom[2], bottom[1], top[0], top[2], top[1]]
    return n


def build_hexahedron(oriented, P):
    """Build a hexahedron connectivity with positive volume orientation."""
    bottom = next(f for f in oriented if len(f) == 4)
    top = match_top(bottom, oriented)
    if top is None:
        return None
    n = list(bottom) + top
    p = [P[i] for i in n]
    if triple(p[1] - p[0], p[3] - p[0], p[4] - p[0]) < 0:
        n = [bottom[0], bottom[3], bottom[2], bottom[1], top[0], top[3], top[2], top[1]]
    return n


def build_polyhedra(poly_cells):
    """Split general polyhedra by unique node count -> polyhedronN CellBlocks."""
    by_n = defaultdict(list)
    for oriented in poly_cells:
        n_nodes = len(set().union(*oriented))
        by_n[n_nodes].append([list(f) for f in oriented])
    cells = []
    for n_nodes, polys in by_n.items():
        data = np.empty(len(polys), dtype=object)
        for i, p in enumerate(polys):
            data[i] = [np.array(f, dtype=int) for f in p]
        cells.append(CellBlock(f"polyhedron{n_nodes}", data))
    return cells


def reconstruct_cell(oriented, P):
    """
    Classify a cell by (n_faces, n_points).

    Returns (meshio_type, connectivity) where:
      - for standard types : connectivity is a flat list of point ids
      - for 'polyhedron'   : connectivity is the list of outward-oriented faces
    """
    n_faces = len(oriented)
    n_pts = len(set().union(*oriented))

    if n_faces == 4 and n_pts == 4:
        return "tetra", build_tetra(oriented, P)
    if n_faces == 5 and n_pts == 5:
        return "pyramid", build_pyramid(oriented, P)
    if n_faces == 5 and n_pts == 6:
        return "wedge", build_wedge(oriented, P)
    if n_faces == 6 and n_pts == 8:
        return "hexahedron", build_hexahedron(oriented, P)

    # General polyhedron: keep outward-oriented faces
    return "polyhedron", oriented


def polygon_from_edges(edges, P):
    """Chain a 2-D cell's boundary edges (``(a, b)`` pairs) into one ring of
    point ids, counter-clockwise in the xy plane. ``None`` when the edges do
    not close a single loop."""
    if len(edges) < 3:
        return None
    nxt = defaultdict(list)
    for a, b in edges:
        nxt[a].append(b)
        nxt[b].append(a)
    if any(len(v) != 2 for v in nxt.values()):
        return None
    start = edges[0][0]
    ring = [start, edges[0][1]]
    while len(ring) < len(nxt):
        a, b = nxt[ring[-1]]
        c = a if a != ring[-2] else b
        if c == start:
            return None
        ring.append(c)
    if len(ring) != len(edges) or ring[0] not in nxt[ring[-1]]:
        return None
    x = [float(P[i][0]) for i in ring]
    y = [float(P[i][1]) for i in ring]
    area = sum(
        x[i] * y[(i + 1) % len(ring)] - x[(i + 1) % len(ring)] * y[i]
        for i in range(len(ring))
    )
    return ring if area >= 0 else [ring[0]] + ring[:0:-1]
