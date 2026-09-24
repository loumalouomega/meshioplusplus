#!/usr/bin/env python3
"""Regenerate the Z88 fixtures under ``tests/python/meshes/z88/``.

Z88OS (GPL-2) is not a dependency and its example decks are not copied: the
structure files are written here from the Z88 manual's layouts, each element's
nodes spelled out by position rather than taken from meshio++'s tables:

* hexahedra (types 1 and 10) number the face the manual draws on top first:
  corners 1-4 counter-clockwise at the top (seen from +z), 5-8 below them; the
  hex20 mid-edge nodes 9-12 run round the top face, 13-16 round the bottom and
  17-20 down the vertical edges 1-5 ... 4-8;
* tet10 (type 16): corners 1-4, then the mid-edge nodes of 1-2, 2-3, 3-1, 2-4,
  3-4, 1-4;
* the plane elements (types 7 and 14): corners counter-clockwise, then the
  mid-edge nodes in the same order.

Four decks, one directory each (Z88's file names are fixed):

* ``cantilever/``: a 100 x 10 x 10 cantilever of five hex20 (type 10), clamped
  at x = 0 and loaded at x = 100, with the ``z88i2``/``z88i5``/``z88mat``/
  ``z88elp``/``z88int``/``z88man`` inputs Z88R needs. Its ``z88o2.txt`` and
  ``z88o3.txt`` are **not** generated here: they are Z88OS V15's own output
  for this deck (``z88r -t -choly``), committed as produced -- see README.md.
* ``tets/``: one tet10 and one tet4 sharing a face.
* ``plate_v13/``: a 2-D plate of one quad8 (type 7) and one tri6 (type 14) with
  the Z88 <= V13 / Aurora V1 header (material count and flags on line 1) and
  its material line after the elements.
* ``polar/``: two tri6 in cylindrical input (``KFLAG = 1``: r, phi in degrees).

    python tools/gen_z88_fixtures.py
"""

import pathlib

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "z88"
)


class Deck:
    def __init__(self, ndim):
        self.ndim = ndim
        self.ids = {}
        self.coords = []
        self.elements = []

    def node(self, p):
        key = tuple(round(float(x), 9) for x in p)
        if key not in self.ids:
            self.ids[key] = len(self.coords) + 1
            self.coords.append(key)
        return self.ids[key]

    def element(self, code, points):
        self.elements.append((code, [self.node(p) for p in points]))

    def header(self, dof):
        return self.ndim, len(self.coords), len(self.elements), dof * len(self.coords)

    def body(self, dof):
        lines = []
        for k, c in enumerate(self.coords, start=1):
            xs = "".join(f"  {v:+.8E}" for v in c[: self.ndim])
            lines.append(f"{k:6d}  {dof}{xs}\n")
        for k, (code, nodes) in enumerate(self.elements, start=1):
            lines.append(f"{k:6d}{code:6d}\n")
            lines.append("".join(f"{n:6d}" for n in nodes) + "\n")
        return "".join(lines)


def hex20(x0, x1, y0, y1, z0, z1):
    """The 20 positions of a type-10 hex in Z88's order (top face first)."""
    top = [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    bot = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)]

    def mid(a, b):
        return tuple((p + q) / 2 for p, q in zip(a, b))

    ring = lambda f: [mid(f[k], f[(k + 1) % 4]) for k in range(4)]  # noqa: E731
    return top + bot + ring(top) + ring(bot) + [mid(top[k], bot[k]) for k in range(4)]


def cantilever():
    d = Deck(3)
    for e in range(5):
        d.element(10, hex20(20.0 * e, 20.0 * (e + 1), 0.0, 10.0, 0.0, 10.0))
    ndim, nn, ne, ndof = d.header(3)
    head = f"{ndim:5d}{nn:6d}{ne:6d}{ndof:7d}{0:6d}   meshio++ Z88 fixture\n"
    out = OUT / "cantilever"
    out.mkdir(parents=True, exist_ok=True)
    (out / "z88i1.txt").write_text(head + d.body(3))
    # Clamp every node at x = 0; share 1000 N in -z over the nodes at x = 100.
    fixed = [k for k, c in enumerate(d.coords, start=1) if c[0] == 0.0]
    loaded = [k for k, c in enumerate(d.coords, start=1) if c[0] == 100.0]
    rows = [f"{k:6d}  {dof}  2  0\n" for k in fixed for dof in (1, 2, 3)]
    rows += [f"{k:6d}  3  1  {-1000.0 / len(loaded):+.6E}\n" for k in loaded]
    (out / "z88i2.txt").write_text(f"{len(rows)}\n" + "".join(rows))
    (out / "z88i5.txt").write_text("0\n")
    (out / "z88mat.txt").write_text(f"1\n1 {ne} 51.txt\n")
    (out / "51.txt").write_text("210000 0.3\n")
    (out / "z88elp.txt").write_text(f"1\n1 {ne} 0 0 0 0 0 0 0\n")
    (out / "z88int.txt").write_text(f"1\n1 {ne} 3 3\n")


def tets():
    d = Deck(3)
    c = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]

    def mid(a, b):
        return tuple((p + q) / 2 for p, q in zip(c[a], c[b]))

    # Z88's tet10: mid-edge nodes of 1-2, 2-3, 3-1, 2-4, 3-4, 1-4.
    d.element(
        16, c + [mid(0, 1), mid(1, 2), mid(2, 0), mid(1, 3), mid(2, 3), mid(0, 3)]
    )
    # A tet4 on the face 2-3-4, apex at (1, 1, 1) outside it.
    d.element(17, [c[1], c[2], c[3], (1.0, 1.0, 1.0)])
    ndim, nn, ne, ndof = d.header(3)
    out = OUT / "tets"
    out.mkdir(parents=True, exist_ok=True)
    (out / "z88i1.txt").write_text(
        f"{ndim:5d}{nn:6d}{ne:6d}{ndof:7d}{0:6d}\n" + d.body(3)
    )


def plate_v13():
    d = Deck(2)
    q = [(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)]
    mq = [
        ((q[k][0] + q[(k + 1) % 4][0]) / 2, (q[k][1] + q[(k + 1) % 4][1]) / 2)
        for k in range(4)
    ]
    d.element(7, q + mq)
    t = [(2.0, 0.0), (3.0, 0.5), (2.0, 1.0)]
    mt = [
        ((t[k][0] + t[(k + 1) % 3][0]) / 2, (t[k][1] + t[(k + 1) % 3][1]) / 2)
        for k in range(3)
    ]
    d.element(14, t + mt)
    ndim, nn, ne, ndof = d.header(2)
    # ndim nnodes nelem ndof nmat kflag ibflag ipflag iqflag
    head = f"{ndim:5d}{nn:6d}{ne:6d}{ndof:7d}{1:6d}{0:6d}{0:6d}{0:6d}{0:6d}\n"
    material = f"{1:6d}{ne:6d}  2.10000E+05  3.00000E-01{3:6d}  1.00000E+00\n"
    out = OUT / "plate_v13"
    out.mkdir(parents=True, exist_ok=True)
    (out / "z88i1.txt").write_text(head + d.body(2) + material)


def polar():
    d = Deck(2)
    # r, phi (degrees): a quarter ring 1 <= r <= 2, 0 <= phi <= 90.
    t1 = [(1.0, 0.0), (2.0, 0.0), (2.0, 90.0), (1.5, 0.0), (2.0, 45.0), (1.5, 45.0)]
    t2 = [(1.0, 0.0), (2.0, 90.0), (1.0, 90.0), (1.5, 45.0), (1.5, 90.0), (1.0, 45.0)]
    d.element(14, t1)
    d.element(14, t2)
    ndim, nn, ne, ndof = d.header(2)
    out = OUT / "polar"
    out.mkdir(parents=True, exist_ok=True)
    (out / "z88i1.txt").write_text(
        f"{ndim:5d}{nn:6d}{ne:6d}{ndof:7d}{1:6d}\n" + d.body(2)
    )


def main():
    cantilever()
    tets()
    plate_v13()
    polar()


if __name__ == "__main__":
    main()
