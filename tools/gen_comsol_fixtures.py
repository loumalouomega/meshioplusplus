#!/usr/bin/env python3
"""Regenerate the COMSOL fixtures under ``tests/python/meshes/comsol/``.

The files are written here, not by meshio++, so the reader is never only checked
against its own writer. Element nodes are listed in COMSOL's own numbering,
built from the rule of COMSOL's "Mesh Element Numbering Conventions" rather than
taken from ``meshioplusplus._node_order``: the corners in tensor order (x
fastest), then every other node of the element's quadratic lattice in
lexicographic (z, y, x) order. The tests check the reader puts every mid-edge
node on its edge's midpoint.

* ``two_domains.mphtxt`` / ``.mphbin`` -- Mesh version 4: a tetrahedron (domain
  1) under a prism (domain 2), their boundary triangles and quads, edges and
  points, and three Selections: a label with blanks, one with a ``#`` (a comment
  character outside strings) and a boundary selection. The ``.mphbin`` is the
  same values in COMSOL's binary serialisation (little-endian int32/float64,
  one int32 per string character).
* ``quadratic.mphtxt`` -- version 4: one each of ``tet2``, ``pyr2``,
  ``prism2``, ``hex2``, ``tri2``, ``quad2`` and ``edg2``.
* ``legacy_v2.mphtxt`` -- Mesh version 2 (lowest vertex index 1): a ``hex2``
  with ``quad2`` boundary elements whose parameter records hold three values
  per node, ``edg2`` elements with one, up/down pairs.
* ``two_objects.mphtxt`` -- two Mesh objects and a Selection of the second.

    python tools/gen_comsol_fixtures.py
"""

import pathlib
import struct

import numpy as np

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "comsol"
)

H = (0.0, 0.5, 1.0)


def lattice(kind):
    """Reference coordinates in COMSOL order for a COMSOL element type."""
    if kind == "edg2":
        verts = [(0, 0, 0), (1, 0, 0)]
        pts = [(x, 0, 0) for x in H]
    elif kind == "tri2":
        verts = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
        pts = [(x, y, 0) for y in H for x in H if x + y <= 1]
    elif kind == "quad2":
        verts = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)]
        pts = [(x, y, 0) for y in H for x in H]
    elif kind == "tet2":
        verts = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
        pts = [(x, y, z) for z in H for y in H for x in H if x + y + z <= 1]
    elif kind == "prism2":
        verts = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)]
        pts = [(x, y, z) for z in H for y in H for x in H if x + y <= 1]
    elif kind == "hex2":
        verts = [(x, y, z) for z in (0, 1) for y in (0, 1) for x in (0, 1)]
        pts = [(x, y, z) for z in H for y in H for x in H]
    elif kind == "pyr2":
        # the base lattice, then the four half-way points to the apex
        apex = (0.5, 0.5, 1.0)
        verts = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), apex]
        pts = [(x, y, 0) for y in H for x in H]
        pts += [tuple((np.array(v) + apex) / 2) for v in verts[:4]]
    else:
        raise ValueError(kind)
    verts = [tuple(float(c) for c in v) for v in verts]
    rest = [tuple(float(c) for c in p) for p in pts]
    return verts + [p for p in rest if p not in verts]


class Values:
    """A COMSOL value stream, rendered as text or as binary."""

    def __init__(self):
        self.items = []  # (kind, value, comment); kind in "i", "d", "s", "#", "\n"

    def i(self, *values, comment=None):
        self.items.append(("i", values, comment))

    def d(self, *values):
        self.items.append(("d", values, None))

    def s(self, text, comment=None):
        self.items.append(("s", text, comment))

    def comment(self, text):
        self.items.append(("#", text, None))

    def text(self):
        out = ["# Created by tools/gen_comsol_fixtures.py (COMSOL layout)\n"]
        for kind, v, comment in self.items:
            if kind == "i":
                line = " ".join(str(x) for x in v)
            elif kind == "d":
                line = " ".join(repr(float(x)) for x in v)
            elif kind == "s":
                line = f"{len(v)} {v}"
            else:
                out.append(f"# {v}\n")
                continue
            out.append(line + (f" # {comment}" if comment else "") + "\n")
        return "".join(out)

    def binary(self):
        out = []
        for kind, v, _ in self.items:
            if kind == "i":
                out.append(struct.pack(f"<{len(v)}i", *v))
            elif kind == "d":
                out.append(struct.pack(f"<{len(v)}d", *v))
            elif kind == "s":
                out.append(struct.pack(f"<{len(v) + 1}i", len(v), *map(ord, v)))
        return b"".join(out)


def header(val, tags):
    val.i(0, 1)
    val.i(len(tags), comment="number of tags")
    for t in tags:
        val.s(t)
    val.i(len(tags), comment="number of types")
    for _ in tags:
        val.s("obj")


def mesh_object(val, points, types, version=4, lowest=0, params=None, updown=None):
    """``types``: [(name, element rows (0-based), entity indices)]."""
    val.comment("--------- Object ----------")
    val.i(0, 0, 1)
    val.s("Mesh", comment="class")
    val.i(version, comment="version")
    val.i(len(points[0]), comment="sdim")
    val.i(len(points), comment="number of mesh vertices")
    val.i(lowest, comment="lowest mesh vertex index")
    val.comment("Mesh vertex coordinates")
    for p in points:
        val.d(*p)
    val.i(len(types), comment="number of element types")
    for name, rows, geom in types:
        val.comment("Type")
        val.s(name, comment="type name")
        val.i(len(rows[0]), comment="number of vertices per element")
        val.i(len(rows), comment="number of elements")
        val.comment("Elements")
        for r in rows:
            val.i(*[x + lowest for x in r])
        if version < 4:
            per, records = (params or {}).get(name, (len(rows[0]), []))
            val.i(per, comment="number of parameter values per element")
            val.i(len(records), comment="number of parameters")
            for rec in records:
                val.d(*rec)
        val.i(len(geom), comment="number of geometric entity indices")
        for g in geom:
            val.i(g)
        if version < 4:
            pairs = (updown or {}).get(name, [])
            val.i(len(pairs), comment="number of up/down pairs")
            for p in pairs:
                val.i(*p)


def selection(val, label, mesh_tag, dim, entities):
    val.comment("--------- Selection ----------")
    val.i(0, 0, 1)
    val.s("Selection", comment="class")
    val.i(0, comment="Version")
    val.s(label, comment="Label")
    val.s(mesh_tag, comment="Geometry/mesh tag")
    val.i(dim, comment="Dimension")
    val.i(len(entities), comment="Number of entities")
    for e in entities:
        val.i(e)


def two_domains():
    # a prism (domain 2) standing on a tetrahedron (domain 1); vertex order
    # scrambled so element order and point order differ
    pts = [
        (0.0, 1.0, 0.0),  # 0
        (0.0, 0.0, 0.0),  # 1
        (1.0, 0.0, 0.0),  # 2
        (0.0, 0.0, -1.0),  # 3  tet apex
        (0.0, 0.0, 1.0),  # 4
        (1.0, 0.0, 1.0),  # 5
        (0.0, 1.0, 1.0),  # 6
    ]
    tet = [[1, 0, 2, 3]]  # positive volume in COMSOL (= VTK) order
    prism = [[1, 2, 0, 4, 5, 6]]
    tris = [[1, 2, 0], [4, 5, 6], [1, 2, 3]]
    quads = [[1, 2, 4, 5]]  # COMSOL quad: tensor order
    edges = [[1, 2], [2, 0]]
    verts = [[1], [3]]
    val = Values()
    header(val, ["mesh1", "mesh1_sel1", "mesh1_sel2", "mesh1_sel3"])
    mesh_object(
        val,
        pts,
        [
            ("vtx", verts, [0, 3]),
            ("edg", edges, [0, 1]),
            ("tri", tris, [0, 2, 1]),
            ("quad", quads, [3]),
            ("tet", tet, [1]),
            ("prism", prism, [2]),
        ],
    )
    selection(val, "Lower Part", "mesh1", 3, [1])
    selection(val, "Part #2", "mesh1", 3, [2])
    selection(val, "caps", "mesh1", 2, [0, 1])
    return val


def quadratic():
    points = []
    types = []
    for k, (name, shift) in enumerate(
        [
            ("tet2", (0, 0, 0)),
            ("pyr2", (3, 0, 0)),
            ("prism2", (6, 0, 0)),
            ("hex2", (9, 0, 0)),
            ("tri2", (0, 3, 0)),
            ("quad2", (3, 3, 0)),
            ("edg2", (6, 3, 0)),
        ]
    ):
        base = len(points)
        points += [tuple(c + s for c, s in zip(p, shift)) for p in lattice(name)]
        types.append((name, [list(range(base, len(points)))], [0 if k >= 4 else 1]))
    val = Values()
    header(val, ["mesh1"])
    mesh_object(val, points, types)
    return val


def legacy_v2():
    lat = lattice("hex2")
    points = lat
    index = {p: k for k, p in enumerate(points)}
    hexa = [list(range(27))]
    # two quad2 faces (z = 0 and z = 1) in their own lattice order
    faces = []
    for z in (0.0, 1.0):
        faces.append([index[(x, y, z)] for x, y, _ in lattice("quad2")])
    edg = [[index[(0.0, 0.0, 0.0)], index[(1.0, 0.0, 0.0)], index[(0.5, 0.0, 0.0)]]]
    vtx = [[index[(0.0, 0.0, 0.0)]]]
    # parameters: three values per node on faces (like COMSOL's 3-D boundary
    # records), one per node on edges
    face_params = [
        [v for p in lattice("quad2") for v in (p[0], p[1], 5.0)] for _ in faces
    ]
    edge_params = [[0.0, 1.0, 0.5]]
    val = Values()
    header(val, ["mesh1"])
    mesh_object(
        val,
        points,
        [
            ("vtx", vtx, [0]),
            ("edg2", edg, [0]),
            ("quad2", faces, [4, 5]),
            ("hex2", hexa, [1]),
        ],
        version=2,
        lowest=1,
        params={"vtx": (1, []), "edg2": (3, edge_params), "quad2": (9, face_params)},
        updown={"quad2": [[0, 1], [0, 1]]},
    )
    return val


def two_objects():
    val = Values()
    header(val, ["mesh1", "mesh2", "mesh2_sel1"])
    mesh_object(val, [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], [("tri", [[0, 1, 2]], [1])])
    mesh_object(
        val,
        [(2.0, 0.0), (3.0, 0.0), (2.0, 1.0), (3.0, 1.0)],
        [("tri", [[0, 1, 2], [1, 3, 2]], [1, 2])],
    )
    selection(val, "right", "mesh2", 2, [2])
    return val


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    val = two_domains()
    (OUT / "two_domains.mphtxt").write_text(val.text())
    (OUT / "two_domains.mphbin").write_bytes(val.binary())
    (OUT / "quadratic.mphtxt").write_text(quadratic().text())
    (OUT / "legacy_v2.mphtxt").write_text(legacy_v2().text())
    (OUT / "two_objects.mphtxt").write_text(two_objects().text())
    print(f"wrote {sorted(p.name for p in OUT.iterdir())}")


if __name__ == "__main__":
    main()
