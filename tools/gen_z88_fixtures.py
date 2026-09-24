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

Eight decks, one directory each (Z88's file names are fixed):

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
* ``frame/``, ``plate/``, ``shell/``, ``torus/``: 3-D beams and a truss, quad8
  plates, tri6 shells and axisymmetric quad8, each with the inputs Z88R needs
  and, committed as produced, Z88R's ``z88o2``/``z88o3``/``z88o4`` output.

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
        """`dof` is one count for every node, or a list of per-node counts."""
        lines = []
        for k, c in enumerate(self.coords, start=1):
            xs = "".join(f"  {v:+.8E}" for v in c[: self.ndim])
            d = dof[k - 1] if isinstance(dof, list) else dof
            lines.append(f"{k:6d}  {d}{xs}\n")
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


def _quad8(corners):
    """Corners counter-clockwise, then the mid-edge nodes in the same order."""
    mids = [
        tuple((a + b) / 2 for a, b in zip(corners[k], corners[(k + 1) % 4]))
        for k in range(4)
    ]
    return list(corners) + mids


def _deck(name, d, dof, header_dof, rows, elp, integration, tail=""):
    """Write `d` and the inputs Z88R needs: `rows` are `z88i2.txt` rows (node,
    dof, flag, value), `elp` and `integration` one row each for all elements."""
    ndim, nn, ne, _ = d.header(1)
    head = (
        f"{ndim:5d}{nn:6d}{ne:6d}{header_dof:7d}{0:6d}   meshio++ Z88 fixture{tail}\n"
    )
    out = OUT / name
    out.mkdir(parents=True, exist_ok=True)
    (out / "z88i1.txt").write_text(head + d.body(dof))
    lines = [f"{n:6d}  {k}  {flag}  {v:+.6E}\n" for n, k, flag, v in rows]
    (out / "z88i2.txt").write_text(f"{len(lines)}\n" + "".join(lines))
    (out / "z88i5.txt").write_text("0\n")
    (out / "z88mat.txt").write_text(f"1\n1 {ne} 51.txt\n")
    (out / "51.txt").write_text("210000 0.3\n")
    (out / "z88elp.txt").write_text(f"1\n1 {ne} {elp}\n")
    (out / "z88int.txt").write_text(f"1\n1 {ne} {integration}\n")


def frame():
    """A 3-D portal: two columns and a girder of beams (type 2, 6 dof) and a
    diagonal truss (type 4); clamped feet, a sway load at the top."""
    d = Deck(3)
    a, b = (0.0, 0.0, 0.0), (100.0, 0.0, 0.0)
    c, e = (0.0, 0.0, 100.0), (100.0, 0.0, 100.0)
    for p, q in ((a, c), (b, e), (c, e)):
        d.element(2, [p, q])
    d.element(4, [a, e])
    rows = [(d.ids[f], k, 2, 0.0) for f in (a, b) for k in range(1, 7)]
    rows.append((d.ids[c], 1, 1, 1000.0))
    # A 10 x 10 section: area, Iyy, e_y, Izz, e_z, It, Wt.
    elp = "100 833.333 5 833.333 5 1406 208"
    _deck("frame", d, 6, 6 * len(d.coords), rows, elp, "0 0")


def plate():
    """Four quad8 Reissner-Mindlin plates (type 20, 3 dof: w and two
    rotations) on [0, 100]^2, clamped at x = 0, a point load at (100, 100)."""
    d = Deck(2)
    for i in range(2):
        for j in range(2):
            x0, y0 = 50.0 * i, 50.0 * j
            q = [(x0, y0), (x0 + 50, y0), (x0 + 50, y0 + 50), (x0, y0 + 50)]
            d.element(20, _quad8(q))
    rows = [
        (k, dof, 2, 0.0)
        for k, c in enumerate(d.coords, start=1)
        if c[0] == 0.0
        for dof in (1, 2, 3)
    ]
    rows.append((d.ids[(100.0, 100.0)], 1, 1, -100.0))
    _deck("plate", d, 3, 3 * len(d.coords), rows, "10 0 0 0 0 0 0", "3 3")


def _tri6(corners):
    """Corners counter-clockwise, then the mid-edge nodes of 1-2, 2-3, 3-1."""
    mids = [
        tuple((a + b) / 2 for a, b in zip(corners[k], corners[(k + 1) % 3]))
        for k in range(3)
    ]
    return list(corners) + mids


def shell():
    """Four tri6 shells (type 24, 6 dof) flat in the plane z = 0 of a 3-D file,
    clamped at x = 0 and pulled and bent at x = 200."""
    d = Deck(3)
    for i in range(2):
        x0 = 100.0 * i
        a, b = (x0, 0.0, 0.0), (x0 + 100, 0.0, 0.0)
        c, e = (x0 + 100, 50.0, 0.0), (x0, 50.0, 0.0)
        d.element(24, _tri6([a, b, c]))
        d.element(24, _tri6([a, c, e]))
    rows = [
        (k, dof, 2, 0.0)
        for k, c in enumerate(d.coords, start=1)
        if c[0] == 0.0
        for dof in range(1, 7)
    ]
    rows += [
        (k, dof, 1, v)
        for k, c in enumerate(d.coords, start=1)
        if c[0] == 200.0
        for dof, v in ((1, 1000.0), (3, -10.0))
    ]
    _deck("shell", d, 6, 6 * len(d.coords), rows, "5 0 0 0 0 0 0", "3 3")


def torus():
    """Two axisymmetric quad8 (type 8) in (r, z), 50 <= r <= 60: held axially
    at z = 0 and pushed outwards at r = 50."""
    d = Deck(2)
    for j in range(2):
        z0 = 10.0 * j
        d.element(8, _quad8([(50.0, z0), (60.0, z0), (60.0, z0 + 10), (50.0, z0 + 10)]))
    rows = [(k, 2, 2, 0.0) for k, c in enumerate(d.coords, start=1) if c[1] == 0.0]
    rows += [(k, 1, 1, 100.0) for k, c in enumerate(d.coords, start=1) if c[0] == 50.0]
    _deck("torus", d, 2, 2 * len(d.coords), rows, "0 0 0 0 0 0 0", "3 3")


def main():
    cantilever()
    frame()
    plate()
    shell()
    torus()
    tets()
    plate_v13()
    polar()


if __name__ == "__main__":
    main()
