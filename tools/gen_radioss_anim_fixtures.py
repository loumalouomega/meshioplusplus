#!/usr/bin/env python3
"""Regenerate the OpenRadioss animation fixtures under
``tests/python/meshes/radioss_anim/``.

OpenRadioss is not run here: the files are written value by value in the layout
its ``anim_to_vtk`` converter reads (big-endian, magic ``0x542C``; OpenRadioss
Tools, MIT), and ``anim_to_vtk`` itself read them outside the repository. One
small model, two states (``cubeA001`` at t = 0, ``cubeA002`` at t = 0.5 with
moved nodes and one element deleted), with every section the reader knows:

* 2-D: a quad shell and a triangle written as a facet with a repeated node, a
  nodal scalar, two nodal vectors, an element scalar and a 2-D tensor;
* 3-D: a hexahedron and a tetrahedron and a wedge written as degenerate bricks,
  an element scalar shared with the 2-D section and a stress tensor;
* 1-D: two beams with an element scalar and a 9-value force/moment set;
* SPH: two particles with a scalar;
* masses, node and element ids, the part/subset/material/property hierarchy
  and time-history lists (which the reader skips).

    python tools/gen_radioss_anim_fixtures.py
"""

import pathlib
import struct

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "radioss_anim"
)


class Out:
    def __init__(self):
        self.b = bytearray()

    def ints(self, *v):
        self.b += struct.pack(f">{len(v)}i", *v)

    def floats(self, *v):
        self.b += struct.pack(f">{len(v)}f", *v)

    def text(self, s, n):
        raw = s.encode()[: n - 1]
        self.b += raw + b"\0" * (n - len(raw))

    def raw(self, data):
        self.b += data


# 15 nodes: a unit hex (0-7) whose top face is the quad shell's; a tetra and a
# wedge beside it on +x; the triangle on the hex's and tetra's top edge; two SPH
# particles on nodes 12 and 13.
BASE = [
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
    (2, 0, 0),
    (2, 1, 0),
    (2, 0, 1),
    (2, 1, 1),
    (3, 0, 0),
    (3, 1, 0),
    (3, 0, 1),
]
N = len(BASE)


def state(time, shift, deleted_brick):
    o = Out()
    o.ints(0x542C)
    o.floats(time)
    for t in ("Time=", "ModAnim", "Radioss Run="):
        o.text(t, 81)
    flags = [1, 1, 1, 1, 1, 1, 0, 1, 0, 0]
    o.ints(*flags)
    # 2-D: two facets (a quad on the hex top, a triangle 4-5-10 written 4 5 10 10)
    facets = [(4, 5, 6, 7), (5, 10, 11, 11)]
    nfun, nefun, nvec, nten, nskew = 1, 1, 2, 1, 0
    o.ints(N, len(facets), 2, nfun, nefun, nvec, nten, nskew)
    for k, (x, y, z) in enumerate(BASE):
        o.floats(x + shift * k / N, y, z)
    for f in facets:
        o.ints(*f)
    o.raw(b"\xff\xff")  # alive (the facets write 0xFF)
    o.ints(1, 2)  # part ends: part "1" has facet 0, part "2" facet 1
    o.text("       11:top shell", 50)
    o.text("       12:side tria", 50)
    o.raw(struct.pack(f">{3 * N}H", *([3000, 0, 0] * N)))  # normals
    o.text("Temperature", 81)
    o.text("Thickness", 81)
    o.floats(*[20.0 + k for k in range(N)])
    o.floats(1.5, 2.5)
    o.text("Velocity", 81)
    o.text("Displacement", 81)
    o.floats(*[v for k in range(N) for v in (k, 0.0, 0.0)])
    o.floats(*[v for k in range(N) for v in (shift * k / N, 0.0, 0.0)])
    o.text("Stress", 81)
    o.floats(10.0, 20.0, 5.0, 11.0, 21.0, 6.0)  # xx yy xy per facet
    o.floats(0.1, 0.2)  # element masses
    o.floats(*[0.01 * (k + 1) for k in range(N)])  # nodal masses
    o.ints(*[100 + k for k in range(N)])  # node ids
    o.ints(501, 502)  # facet ids
    o.ints(1, 1)  # subsets
    o.ints(7, 7)  # materials
    o.ints(3, 4)  # properties
    # 3-D: a hex, a tetra (1 2 3 3 4 4 4 4 pattern) and a wedge (1 2 3 3 5 6 7 7)
    bricks = [
        (0, 1, 2, 3, 4, 5, 6, 7),
        (1, 8, 9, 9, 10, 10, 10, 10),
        (8, 12, 9, 9, 10, 14, 11, 11),
    ]
    o.ints(len(bricks), 2, 1, 1)
    for b in bricks:
        o.ints(*b)
    o.raw(bytes([0 if k == deleted_brick else 1 for k in range(len(bricks))]))
    o.ints(1, 3)  # part "21" has brick 0, part "22" bricks 1-2
    o.text("       21:block", 50)
    o.text("       22:wedges", 50)
    o.text("Thickness", 81)  # the same name as the 2-D scalar
    o.floats(7.0, 8.0, 9.0)
    o.text("Stress", 81)  # the same name as the 2-D tensor
    o.floats(*[float(10 * b + c) for b in range(3) for c in range(6)])
    o.floats(1.0, 0.5, 0.5)  # masses
    o.ints(1001, 1002, 1003)  # ids
    o.ints(1, 1)
    o.ints(8, 9)
    o.ints(5, 6)
    # 1-D: two beams with a scalar and a force/moment set, skews
    beams = [(0, 8), (8, 12)]
    o.ints(len(beams), 1, 1, 1, 1)
    for b in beams:
        o.ints(*b)
    o.raw(b"\x01\x01")
    o.ints(2)
    o.text("       31:frame", 50)
    o.text("Plastic Strain", 81)
    o.floats(0.01, 0.02)
    o.text("Forces-Moments", 81)
    o.floats(*[float(b * 9 + c) for b in range(2) for c in range(9)])
    o.ints(0, 0)  # skews
    o.floats(0.3, 0.4)
    o.ints(2001, 2002)
    o.ints(1)
    o.ints(7)
    o.ints(2)
    # hierarchy: one subset holding every part, two materials, three properties
    o.ints(1)
    o.text("model", 50)
    o.ints(0, 0)  # parent, no child subsets
    o.ints(2, 0, 1)  # its 2-D parts
    o.ints(2, 0, 1)  # 3-D parts
    o.ints(1, 0)  # 1-D parts
    o.ints(2, 3)
    o.text("steel", 50)
    o.text("rubber", 50)
    o.ints(2, 42)
    for p in ("shell", "solid", "beam"):
        o.text(p, 50)
    o.ints(1, 14, 3)
    # time-history lists: 1 node, 1 facet, 0 bricks, 1 beam
    o.ints(1, 1, 0, 1)
    o.ints(3)
    o.text("tip", 50)
    o.ints(0)
    o.text("shell", 50)
    o.ints(1)
    o.text("beam", 50)
    # SPH: two particles on nodes 12 and 13
    o.ints(2, 1, 1, 0)
    o.ints(12, 13)
    o.raw(b"\x01\x01")
    o.ints(2)
    o.text("       41:water", 50)
    o.text("Diameter", 81)
    o.floats(0.25, 0.25)
    o.floats(0.05, 0.05)
    o.ints(3001, 3002)
    o.ints(1)
    o.ints(9)
    o.ints(7)
    return bytes(o.b)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "cubeA001").write_bytes(state(0.0, 0.0, -1))
    (OUT / "cubeA002").write_bytes(state(0.5, 0.25, 1))


if __name__ == "__main__":
    main()
