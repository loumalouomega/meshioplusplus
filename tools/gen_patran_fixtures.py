#!/usr/bin/env python3
"""Regenerate the Patran 2 neutral fixtures under ``tests/python/meshes/patran/``.

The cards are written here, not by meshio++, so the reader is never only checked
against its own writer. They follow the packet layouts of the Patran "Interface
to PATRAN 2 Neutral File" guide -- a ``(I2,8I8)`` header card ``IT, ID, IV, KC,
N1..N5`` followed by ``KC`` data cards -- and each quadratic element lists its
mid-edge nodes in the order of the Patran Reference Manual's Element Library,
spelled out below as edge lists rather than taken from
``meshioplusplus._node_order``:

* bar3 ``(1,2)``; tri6 ``(1,2) (2,3) (3,1)``; quad8 ``(1,2) (2,3) (3,4) (4,1)``;
* tet10 ``(1,2) (2,3) (3,1) (1,4) (2,4) (3,4)``;
* pyramid13: the base ring ``(1,2) (2,3) (3,4) (4,1)``, then ``(1,5) ... (4,5)``;
* wedge15: bottom ring ``(1,2) (2,3) (3,1)``, the *vertical* edges ``(1,4) (2,5)
  (3,6)``, then the top ring ``(4,5) (5,6) (6,4)``;
* hex20: bottom ring, the vertical edges ``(1,5) ... (4,8)``, then the top ring.

Every mid-edge node sits at its edge's midpoint, so a test can check the read
cells against meshio++'s own edge tables.

Two files:

* ``mixed_linear.pat``: one hex, wedge, tet and pyramid, two quads and a bar,
  with gaps in the node and element ids, a title and summary, packets the reader
  skips (04 property, 06 load, 10 a nodal temperature), and three components:
  ``SOLIDS`` (the four solids), ``FIXED`` (four nodes) and ``MIXED`` (a quad
  and two nodes). The quads and bar are in no component, so they fall back to
  ``property_<pid>`` regions.
* ``quadratic.pat``: one element of each quadratic shape, each in its own
  property, with no components.

    python tools/gen_patran_fixtures.py
"""

import pathlib

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "patran"
)

# Patran shape code (IV) and the component entity code of each shape.
SHAPE = {"bar": (2, 6), "tri": (3, 7), "quad": (4, 8), "tet": (5, 9)}
SHAPE.update({"pyr": (6, 10), "wedge": (7, 11), "hex": (8, 12)})

# Mid-edge node order of each quadratic shape (1-based corner pairs), from the
# Patran Element Library.
EDGES = {
    "bar": [(1, 2)],
    "tri": [(1, 2), (2, 3), (3, 1)],
    "quad": [(1, 2), (2, 3), (3, 4), (4, 1)],
    "tet": [(1, 2), (2, 3), (3, 1), (1, 4), (2, 4), (3, 4)],
    "pyr": [(1, 2), (2, 3), (3, 4), (4, 1), (1, 5), (2, 5), (3, 5), (4, 5)],
    "wedge": [(1, 2), (2, 3), (3, 1), (1, 4), (2, 5), (3, 6), (4, 5), (5, 6), (6, 4)],
    "hex": [(1, 2), (2, 3), (3, 4), (4, 1), (1, 5), (2, 6), (3, 7), (4, 8)]
    + [(5, 6), (6, 7), (7, 8), (8, 5)],
}

# Reference corners (Patran and VTK agree on corner order).
CORNERS = {
    "bar": [(0, 0, 0), (1, 0, 0)],
    "tri": [(0, 0, 0), (1, 0, 0), (0, 1, 0)],
    "quad": [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)],
    "tet": [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
    "pyr": [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1)],
    "wedge": [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)],
    "hex": [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
    + [(0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)],
}


def header(it, ident, iv, kc, *n):
    n = list(n) + [0] * (5 - len(n))
    return f"{it:2d}{ident:8d}{iv:8d}{kc:8d}" + "".join(f"{v:8d}" for v in n) + "\n"


def e16(v):
    return f"{v:16.9E}"


def int_cards(values):
    return "".join(
        "".join(f"{v:8d}" for v in values[k : k + 10]) + "\n"
        for k in range(0, len(values), 10)
    )


class Model:
    def __init__(self, node_step, elem_step):
        self.nodes = []  # (id, xyz)
        self.elements = []  # (id, shape, pid, node ids)
        self.node_step = node_step
        self.elem_step = elem_step

    def node(self, xyz):
        ident = (len(self.nodes) + 1) * self.node_step
        self.nodes.append((ident, xyz))
        return ident

    def element(self, shape, pid, ids):
        ident = (len(self.elements) + 1) * self.elem_step + 100
        self.elements.append((ident, shape, pid, ids))
        return ident

    def solid(self, shape, pid, offset, quadratic=False):
        corners = [tuple(c + o for c, o in zip(p, offset)) for p in CORNERS[shape]]
        ids = [self.node(p) for p in corners]
        if quadratic:
            for a, b in EDGES[shape]:
                pa, pb = corners[a - 1], corners[b - 1]
                ids.append(self.node(tuple((x + y) / 2 for x, y in zip(pa, pb))))
        return self.element(shape, pid, ids)

    def text(self, title, extra="", components=()):
        out = [header(25, 0, 0, 1), title.ljust(80)[:80].rstrip() + "\n"]
        out.append(header(26, 0, 0, 1, len(self.nodes), len(self.elements), 1, 2, 0))
        out.append("23-Sep-26   12:00:00         3.0\n")
        for ident, xyz in self.nodes:
            out.append(header(1, ident, 0, 2))
            out.append("".join(e16(v) for v in xyz) + "\n")
            out.append("1G       6       0       0  000000\n")
        for ident, shape, pid, ids in self.elements:
            kc = 1 + (len(ids) + 9) // 10
            out.append(header(2, ident, SHAPE[shape][0], kc))
            out.append(f"{len(ids):8d}{0:8d}{pid:8d}{0:8d}" + e16(0.0) * 3 + "\n")
            out.append(int_cards(ids))
        out.append(extra)
        for number, name, pairs in components:
            flat = [v for pair in pairs for v in pair]
            out.append(header(21, number, len(flat), 1 + (len(flat) + 9) // 10))
            out.append(name + "\n")
            out.append(int_cards(flat))
        out.append(header(99, 0, 0, 1))
        return "".join(out)


def mixed_linear():
    m = Model(node_step=10, elem_step=7)
    hexa = m.solid("hex", 1, (0, 0, 0))
    wedge = m.solid("wedge", 1, (2, 0, 0))
    tet = m.solid("tet", 1, (4, 0, 0))
    pyr = m.solid("pyr", 1, (6, 0, 0))
    q1 = m.solid("quad", 2, (0, 3, 0))
    q2 = m.solid("quad", 2, (1, 3, 0))
    bar = m.solid("bar", 3, (0, 5, 0))
    first = m.nodes[0][0]
    # Packets the reader skips: a material-less property (04), a distributed
    # load on the hex (06) and a nodal temperature (10).
    extra = (
        header(4, 1, 5, 1)
        + "PSOLID      \n"
        + header(6, hexa, 1, 2, 1, 0)
        + "   0\n"
        + e16(1.0)
        + "\n"
        + header(10, first, 1, 1)
        + e16(300.0)
        + "\n"
    )
    components = [
        (1, "SOLIDS", [(12, hexa), (11, wedge), (9, tet), (10, pyr)]),
        (4, "FIXED", [(5, m.nodes[k][0]) for k in range(4)]),
        (9, "MIXED", [(8, q1), (5, m.nodes[0][0]), (5, m.nodes[1][0])]),
    ]
    assert q2 and bar
    return m.text("mixed linear fixture for meshio++", extra, components)


def quadratic():
    m = Model(node_step=1, elem_step=1)
    for pid, shape in enumerate(
        ["hex", "wedge", "tet", "pyr", "quad", "tri", "bar"], 1
    ):
        m.solid(shape, pid, (3 * pid, 0, 0), quadratic=True)
    return m.text("quadratic fixture for meshio++")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "mixed_linear.pat").write_text(mixed_linear(), newline="\n")
    (OUT / "quadratic.pat").write_text(quadratic(), newline="\n")


if __name__ == "__main__":
    main()
